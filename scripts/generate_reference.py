"""Generate deterministic pymsis(version=0) reference products; no network needed."""
from datetime import datetime, timedelta
from pathlib import Path
import importlib.metadata
import json
import h5py
import numpy as np
import pymsis

ROOT = Path(__file__).resolve().parents[1]
SEED = 82026


def main():
    rng = np.random.default_rng(SEED)
    cases = []
    # Explicit boundaries plus seeded global samples. No downloaded space weather.
    for mode in ("daily", "history"):
        for altitude in (0, 10, 32.5, 62.5, 72.5, 80, 100, 120, 200, 400, 800, 1500, 2000):
            for latitude in (-90, -69, 0, 69, 90):
                index = len(cases)
                date = datetime(2024, 1, 1) + timedelta(days=(index * 37) % 366,
                                                       seconds=(index * 7919) % 86400)
                cases.append(make_case(date, latitude, float((index*53)%360)-180,
                                       altitude, 80 if index%2 else 200,
                                       100 if index%2 else 170, mode, index))
        for index in range(150):
            date = datetime(2023, 1, 1) + timedelta(seconds=int(rng.integers(0, 365*86400)))
            cases.append(make_case(date, float(rng.uniform(-89,89)), float(rng.uniform(-180,180)),
                                   float(rng.uniform(0,2000)), float(rng.uniform(65,250)),
                                   float(rng.uniform(65,250)), mode, index))
    for mode in ("daily", "history"):
        subset = [c for c in cases if c["mode"] == mode]
        output = pymsis.calculate(
            np.array([np.datetime64(c["date"]) for c in subset]),
            np.array([c["longitude_deg"] for c in subset]),
            np.array([c["latitude_deg"] for c in subset]),
            np.array([c["altitude_km"] for c in subset]),
            f107s=np.array([c["f107"] for c in subset]),
            f107as=np.array([c["f107_average"] for c in subset]),
            aps=np.array([c["ap"] for c in subset]),
            version=0, geomagnetic_activity=1 if mode=="daily" else -1,
        ).reshape(-1,11)
        for c, values in zip(subset, output):
            # pymsis order: rho, N2, O2, O, He, H, Ar, N, anomalous O, NO, T.
            c["expected"] = [float(v) if np.isfinite(v) else None for v in values]
    ROOT.joinpath("tests/reference.json").write_text(json.dumps(cases, indent=2)+"\n")
    with h5py.File(ROOT / "tests/reference.h5", "w") as f:
        f.attrs.update(generating_script="scripts/generate_reference.py", seed=SEED,
                       pymsis_version=importlib.metadata.version("pymsis"), model="NRLMSISE-00",
                       units="SI; altitude km; latitude/longitude degrees; temperature K")
        f.create_dataset("cases_json", data=json.dumps(cases))
        f.create_dataset("expected", data=np.array([[np.nan if x is None else x for x in c["expected"]] for c in cases]))
    print(f"Wrote {len(cases)} reference cases against pymsis {importlib.metadata.version('pymsis')}")
    # The linked, unmodified MSIS 2.1 backend uses the same cases, all species.
    for mode in ("daily", "history"):
        subset = [c for c in cases if c["mode"] == mode]
        output = pymsis.calculate(
            np.array([np.datetime64(c["date"]) for c in subset]),
            np.array([c["longitude_deg"] for c in subset]),
            np.array([c["latitude_deg"] for c in subset]),
            np.array([c["altitude_km"] for c in subset]),
            f107s=np.array([c["f107"] for c in subset]),
            f107as=np.array([c["f107_average"] for c in subset]),
            aps=np.array([c["ap"] for c in subset]),version=2.1,
            geomagnetic_activity=1 if mode=="daily" else -1).reshape(-1,11)
        for c, values in zip(subset,output):
            c["expected"]=[float(v) if np.isfinite(v) else None for v in values]
    ROOT.joinpath("tests/reference21.json").write_text(json.dumps(cases,indent=2)+"\n")
    with h5py.File(ROOT/"tests/reference21.h5","w") as f:
        f.attrs.update(generating_script="scripts/generate_reference.py",seed=SEED,
                       pymsis_version=importlib.metadata.version("pymsis"),model="NRLMSIS 2.1")
        f.create_dataset("cases_json",data=json.dumps(cases))
        f.create_dataset("expected",data=np.array([[np.nan if x is None else x for x in c["expected"]] for c in cases]))
    print(f"Wrote {len(cases)} MSIS 2.1 reference cases")


def make_case(date, latitude, longitude, altitude, f107_average, f107, mode, index):
    ap = ([4]*7 if index%3==0 else [25,80,60,40,20,15,8] if index%3==1 else [100,200,150,100,80,60,40])
    return dict(date=date.isoformat(), day_of_year=date.timetuple().tm_yday,
                seconds_of_day=date.hour*3600+date.minute*60+date.second,
                latitude_deg=latitude, longitude_deg=longitude, altitude_km=altitude,
                f107_average=f107_average, f107=f107, ap=ap, mode=mode)


if __name__ == "__main__":
    main()
