# eqlat

**Equivalent latitude (φ_e) from potential vorticity on isentropic surfaces.**

`eqlat` computes equivalent latitude from 2-D fields of Ertel potential vorticity (PV)
with two interchangeable methods. It also includes a small pipeline for reanalysis data
(ERA5, MERRA-2): download the data, interpolate it from pressure levels to isentropic
surfaces, compute φ_e for every time step and level, and merge the results.

- **Piecewise-constant method**: the traditional grid-cell sum.
- **Region of Interest (ROI) method**: contour-based area integration following
  [Añel et al. (2013)](https://doi.org/10.1371/journal.pone.0072970).

Both methods work globally (one formula for both hemispheres). They accept ERA5- and
MERRA-2-style grids as they come and handle missing values, for example where an
isentrope intersects the ground.

---

## Contents

1. [Background](#background)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Batch processing of NetCDF files](#batch-processing-of-netcdf-files)
5. [Reanalysis pipeline (command line)](#reanalysis-pipeline-command-line)
6. [API reference](#api-reference)
7. [Package layout](#package-layout)
8. [Testing](#testing)
9. [Known limitations](#known-limitations)
10. [Reference](#reference)

---

## Background

On an isentropic surface, PV is approximately materially conserved, so PV contours act
like material lines. Equivalent latitude turns a PV value into a latitude: φ_e(q₀) is
the latitude of the circle around the Pole that encloses the **same area** as the
region where PV is greater than q₀.

```
A(q₀) = area where PV > q₀
φ_e   = arcsin(1 − A / (2π R²))          R = 6.371·10⁶ m
```

This convention gives

| A(q₀) | φ_e |
|---|---|
| 0 | 90° N |
| 2πR² (half the sphere) | 0° |
| 4πR² (whole sphere) | 90° S |

so one formula covers both hemispheres: NH PV is positive and SH PV is negative.
The result is a monotonic lookup table PV → φ_e. Applied to every grid point, this
table gives a φ_e field in which transport and mixing barriers (polar vortex edge,
subtropical jet) show up much more clearly than in geographic latitude.

### Piecewise-constant method

Each grid cell is treated as having constant PV over its area. The enclosed area is the
sum of the areas of all cells with PV > q₀:

```
A(q₀) = Σ  a_ij   for all cells with PV_ij > q₀
a_ij  = R² · Δλ · |Δ(sin φ)|
```

The implementation sorts the PV values once and uses a cumulative area sum with
`searchsorted`. This is fast, but the area changes in steps, so the result depends on
the grid resolution.

### ROI method (Añel et al., 2013)

The ROI method integrates the area enclosed by the actual PV contour lines. The
contours are calculated with [`contourpy`](https://contourpy.readthedocs.io), and the
enclosed area is computed with the signed spherical polygon formula

```
A = −R² Σ ½ (sin φ_k + sin φ_{k+1}) · Δλ_k
```

All three contour configurations described in the paper are handled:

| Case | Geometry | Handling |
|---|---|---|
| 1 | Simple closed contour(s) | Area of each closed polygon |
| 2 | Contour split by the dateline | Polygons on both sides are summed |
| 3 | Contour reaching the grid's northern/southern boundary | The signed area formula gives the correct polar cap |

Other details:

- **Hemispheres:** for q₀ ≥ 0, the region `PV ∈ [q₀, PV_max]` is filled (north polar
  cap). For q₀ < 0, the region `PV ∈ [PV_min, q₀]` is filled, and
  `A = A_valid − A_south`.
- **Holes** (low-PV air enclosed by higher PV, e.g. filaments) are handled through the
  winding order of the sub-paths.
- **Missing data:** NaN cells are masked out of the contouring and of the total area.
- **Fallback:** if contouring fails for a threshold, that threshold falls back to the
  piecewise estimate.
- **Grid normalisation:** latitude is flipped to ascending order and longitude is
  shifted to `[0, 360)` internally. Both ERA5 (90 → −90, −180…180) and MERRA-2
  (−90 → 90, −180…179.375) work without preprocessing.

---

## Installation

The package uses a `src/` layout with optional dependency groups:

```bash
git clone <repo-url> eqlat
cd eqlat

pip install -e .              # core algorithms only: numpy, contourpy
pip install -e ".[io]"        # + xarray / NetCDF batch processing, dask, tqdm
pip install -e ".[download]"  # + ERA5 / MERRA-2 / MLS / ozonesonde downloaders
pip install -e ".[all]"       # io + download
pip install -e ".[dev]"       # io + pytest
```

Python ≥ 3.10 is required.

### In a virtual environment

Either with `venv`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[all]"
```

or with conda:

```bash
conda create -n eqlat python=3.12
conda activate eqlat
pip install -e ".[all]"
```

The command-line tools below are available whenever the environment is active.

### Credentials for downloads

| Source | Setup |
|---|---|
| **ERA5** (Copernicus CDS) | Register at <https://cds.climate.copernicus.eu>, then create `~/.cdsapirc`:<br>`url: https://cds.climate.copernicus.eu/api`<br>`key: <your-api-key>` |
| **MERRA-2 / MLS** (NASA Earthdata) | Register at <https://urs.earthdata.nasa.gov>, [authorise GES DISC](https://disc.gsfc.nasa.gov/earthdata-login), then create `~/.netrc`:<br>`machine urs.earthdata.nasa.gov login <user> password <pw>`<br>and run `chmod 600 ~/.netrc` |
| **Ozonesondes** (NOAA GML) | No account needed |

---

## Quick start

### A single 2-D field

```python
import numpy as np
from eqlat import equivalent_latitude_roi, equivalent_latitude_piecewise

# pv  : (nlat, nlon) PV on one isentropic surface, e.g. in PVU (NaN allowed)
# lat : (nlat,)  degrees, ascending or descending
# lon : (nlon,)  degrees, [-180, 180] or [0, 360]

res = equivalent_latitude_roi(pv, lat, lon, n_thresholds=100)

res["pv_thresholds"]   # (n,) PV values q₀
res["eqlat"]           # (n,) equivalent latitude φ_e(q₀) in degrees
res["area"]            # (n,) enclosed area A(q₀) in m²
```

Specific thresholds can be passed as well:

```python
res = equivalent_latitude_piecewise(pv, lat, lon, pv_thresholds=np.linspace(-50, 50, 201))
```

### An equivalent-latitude map

To get φ_e at every grid point, use the `eqlat_field_*` helpers. They build the
PV → φ_e lookup table and interpolate it to each grid point:

```python
from eqlat.roi_fast import eqlat_field_roi
from eqlat.piecewise import eqlat_field_piecewise

eqlat_map = eqlat_field_roi(pv, lat, lon, n_thresholds=200)   # shape (nlat, nlon)
```

### Sanity check

For a zonally symmetric PV field that increases monotonically with latitude, φ_e must
equal geographic latitude:

```python
lat = np.linspace(-89.5, 89.5, 180)
lon = np.arange(0, 360, 1.0)
pv  = np.repeat(np.sin(np.deg2rad(lat))[:, None], lon.size, axis=1)

res = equivalent_latitude_roi(pv, lat, lon, pv_thresholds=[0.5])
res["eqlat"]   # ≈ 30°
```

---

## Batch processing of NetCDF files

`eqlat.batch` (requires the `[io]` extra) processes whole files over all time steps and
levels and returns an `xarray.Dataset`. There are two workflows, depending on the
vertical coordinate of the input.

### 1. PV already on isentropic levels: `process_netcdf`

```python
from eqlat import batch

ds = batch.process_netcdf(
    "era5_pv_theta_20230115.nc",
    method="roi",            # or "piecewise"
    n_thresholds=200,
    n_workers=4,             # parallel processes; 1 = sequential
    theta_slice=slice(0, 5), # optional subset
)
ds["eqlat"]                  # (time, theta, lat, lon), degrees_north
```

### 2. PV and temperature on pressure levels: `process_pressure_netcdf`

This is the standard case for ERA5 pressure-level downloads and MERRA-2 `M2I3NPASM`.
For each column, the function computes θ = T (p₀/p)^κ, interpolates PV linearly in θ
to each requested isentrope, and then computes φ_e.

```python
theta_levels = [320, 340, 350, 360, 380, 400, 450, 500, 550, 600, 700, 800]

ds = batch.process_pressure_netcdf(
    "era5_pv_pressure_20230115.nc",
    theta_levels,
    method="roi",
    include_o3=True,         # also interpolate ozone if the file contains it
    n_workers=4,
)
ds.to_netcdf("era5_eqlat_roi_20230115.nc")
```

The output contains

| Variable | Dims | Description |
|---|---|---|
| `eqlat` | (time, theta, lat, lon) | Equivalent latitude [degrees_north]; attrs record method, `n_thresholds`, source file and variables |
| `pv_isentropic` | (time, theta, lat, lon) | PV interpolated to θ [PVU] |
| `o3_isentropic` | (time, theta, lat, lon) | Ozone interpolated to θ (only if present in the input) |

### Automatic name detection

Variable and dimension names are detected from common conventions. If detection fails,
the error message lists every name that was tried. Each name can be set explicitly
(`pv_var=`, `t_var=`, `o3_var=`, `time_dim=`, `plev_dim=`, `theta_dim=`, `lat_dim=`,
`lon_dim=`).

| Role | Names tried |
|---|---|
| time | `time`, `valid_time` |
| pressure level | `level`, `lev`, `plev`, `pressure`, `pressure_level`, `isobaric`, `hybrid` |
| theta level | `level`, `theta`, `lev`, `plev`, `isentropic_level`, `potential_temperature` |
| latitude | `latitude`, `lat`, `nlat`, `y` |
| longitude | `longitude`, `lon`, `nlon`, `x` |
| PV | `pv`, `epv`, `EPV`, `PV`, `ertel_pv`, `q` |
| temperature | `T`, `t`, `temperature`, `temp`, `ta`, `air_temperature`, `TEMP`, `Temperature` |
| ozone | `o3`, `O3`, `ozone`, `ozone_mass_mixing_ratio`, `mass_fraction_of_ozone_in_air` |

### Error handling

- Each 2-D slice runs in its own `try/except`. A failing slice emits a `RuntimeWarning`
  and stays `NaN`, so a long batch run is not aborted by one bad time step.
- Interpolation and the φ_e computation are guarded separately.
- `xarray` and `tqdm` are optional imports. A clear `ImportError` is raised only when
  they are actually needed. Without `tqdm`, a simple text progress indicator is used.

The shortcuts `process_netcdf_roi`, `process_netcdf_piecewise`,
`process_pressure_netcdf_roi` and `process_pressure_netcdf_piecewise` set `method` for
you.

---

## Reanalysis pipeline (command line)

After installation, the following commands are available. Each one accepts `--help`
and can also be run as `python -m eqlat.<subpackage>.<module>`.

| Command | Purpose |
|---|---|
| `eqlat-download-era5` | ERA5 PV + T on pressure levels (CDS) |
| `eqlat-download-merra2` | MERRA-2 EPV + T on pressure levels (GES DISC, OPeNDAP) |
| `eqlat-download-mls` | Aura MLS ozone + temperature profiles (Earthdata) |
| `eqlat-download-ozonesondes` | NOAA GML ozonesonde profiles → CSV |
| `eqlat-process-fields` | Compute φ_e for a directory of downloaded files |
| `eqlat-merge-reanalysis` | Merge daily output files into one NetCDF |

If you don't pass `--outdir`, the download tools write to subfolders of the directory
in the environment variable `$DATA`, or of the current directory if `$DATA` is not set.

### Example workflow

The example below runs the pipeline for a **single day, 15 January 2023, 12 UTC** on a
local machine. All data ends up in `./data` inside the folder you run it from.

```bash
export DATA=./data

# 0. The date to process (YYYY-MM-DD)
echo 2023-01-15 > dates.txt

# 1. Download ERA5 PV + temperature
eqlat-download-era5 --dates-file dates.txt --hours 12:00 --outdir $DATA/era5
#   → data/era5/era5_pv_pressure_20230115.nc

# 2. Equivalent latitude on 27 isentropes (300–1400 K), piecewise and ROI
eqlat-process-fields --indir $DATA/era5 --outdir $DATA/eqlat --model ERA5 --roi
#   → data/eqlat/era5_eqlat_piecewise_20230115.nc
#   → data/eqlat/era5_eqlat_roi_20230115.nc
```

What each step does:

1. **Download:** fetches PV and temperature on 19 pressure levels (1–350 hPa, 0.25°
   grid) from the Copernicus CDS. CDS requests are queued server-side, so this step can
   take a few minutes even for one day.
2. **Equivalent latitude:** computes θ from T and p, interpolates PV to the isentropes,
   and computes φ_e for every level with both methods.

Look at the result in Python (plotting needs `pip install matplotlib`):

```python
import xarray as xr

ds = xr.open_dataset("data/eqlat/era5_eqlat_roi_20230115.nc")
ds["eqlat"].sel(theta=500).squeeze().plot()          # φ_e map on the 500 K surface
ds["pv_isentropic"].sel(theta=500).squeeze().plot()  # the underlying PV
```

### More than one day

Put one date per line into the dates file, or pass a year instead:

```bash
printf '2023-01-%02d\n' {1..31} > dates.txt                       # all of January
eqlat-download-era5 --dates-file dates.txt --hours 12:00 --outdir $DATA/era5
eqlat-download-era5 2023 --hours 12:00 --outdir $DATA/era5         # a whole year
```

`eqlat-process-fields` then processes every file in `--indir`, one after another.
Existing output files are skipped, so you can stop with `Ctrl+C` and continue later.
Afterwards, join the daily files along the time axis into one file per method:

```bash
eqlat-merge-reanalysis --indir $DATA/eqlat --wildcard_name 'era5_eqlat_roi_*.nc' \
                       --outdir $DATA --outname era5_eqlat_roi_2023.nc
```

Before starting a long run:

- **Disk space.** At 0.25° resolution, one time step of PV + T on 19 levels is roughly
  150 MB uncompressed. A full year at one time step per day is in the range of tens of
  GB, and every additional hour in `--hours` multiplies that. Check the size of the
  single-day file first.
- **Runtime.** The ROI method is considerably slower than the piecewise method. Time
  the single-day example to estimate the full period. If you only need one method,
  leave out `--roi`.

### `eqlat-download-era5`

```
eqlat-download-era5 [YEAR] [--dates-file FILE] [--outdir DIR] [--hours HH:MM ...]
```

- Downloads `potential_vorticity` and `temperature` from
  `reanalysis-era5-pressure-levels` on a 0.25° grid and 19 levels (1–350 hPa).
- Downloads either every day of `YEAR`, or only the dates in `--dates-file` (CSV with a
  `date` column, or a text file with one `YYYY-MM-DD` per line).
- Default output: `$DATA/ERA5_12UTC`. Default hours: all 24, so pass `--hours` to
  keep the download small.
- Files that already exist are skipped.

### `eqlat-download-merra2`

```
eqlat-download-merra2 [YEAR] [--dates-file FILE] [--outdir DIR]
                      [--collection {M2I3NPASM,M2I6NPANA}]
                      [--variables EPV T ...] [--hours 0 6 12 18]
```

- `M2I3NPASM`: 3-hourly instantaneous, 42 pressure levels, 0.5° × 0.625°.
  `M2I6NPANA`: 6-hourly analysis.
- Uses OPeNDAP subsetting, so only the requested variables and time steps are
  downloaded instead of the full daily file of about 1 GB.
- Output: `$DATA/MERRA2/merra2_<VARS>_YYYYMMDD_t<idx>.nc4`.
- The MERRA-2 stream number (100/200/300/400) is chosen from the date.

### `eqlat-download-mls`

```
eqlat-download-mls DATE1 DATE2 [--outdir DIR]
```

Searches and downloads Aura MLS Level 2 `ML2O3` and `ML2T` (v005) granules for the date
range via `earthaccess`. Default output: `$DATA/satellite`.

### `eqlat-download-ozonesondes`

```
eqlat-download-ozonesondes YEAR --url STATION_URL [--outdir DIR]
```

Scrapes a NOAA GML station directory for `.dat` profiles of the given year. Each
profile is parsed (flight number, location, launch coordinates, GMT timestamps; `99999`
→ NaN) and saved as `ozonesonde_<flight>_<location>_<YYYYMMDD>.csv`. Example station
URL:

```
https://gml.noaa.gov/aftp/data/ozwv/Ozonesonde/Boulder,%20Colorado/Native%20Resolution%20(60s,%207s,%201s)/
```

### `eqlat-process-fields`

```
eqlat-process-fields --indir DIR --outdir DIR --model {ERA5,MERRA2} [--roi] [--index N]
```

- Input files are found by pattern: `*pv*` for ERA5, `*EPV*` for MERRA-2.
- The piecewise method always runs. `--roi` also runs the ROI method.
- θ levels: 300–400 K every 10 K, then 425–500 K every 25 K, 550–800 K every 50 K and
  900–1400 K every 100 K (27 levels).
- Output files that already exist are skipped, so interrupted runs can be restarted.
- `--index N` processes only the N-th file (0-based, in sorted order). This is handy
  for a quick test run or to time a single day:

```bash
time eqlat-process-fields --indir $DATA/era5 --outdir $DATA/eqlat --model ERA5 --roi --index 0
```

- Files are processed sequentially on one CPU core. To use several cores, call
  `eqlat.batch.process_pressure_netcdf(..., n_workers=4)` from Python instead (see
  [Batch processing](#batch-processing-of-netcdf-files)).

### `eqlat-merge-reanalysis`

```
eqlat-merge-reanalysis --indir DIR --wildcard_name PATTERN --outdir DIR --outname FILE
```

Opens all matching files lazily with `xr.open_mfdataset(combine="by_coords")`, using
dask with `time` chunks of 10, and writes them out as one file with the `h5netcdf`
engine.

---

## API reference

### `eqlat` (top level)

| Function | Description |
|---|---|
| `equivalent_latitude_piecewise(pv, lat, lon, pv_thresholds=None, n_thresholds=100)` | PV → φ_e table, piecewise-constant method. Returns dict `pv_thresholds`, `eqlat`, `area` |
| `equivalent_latitude_roi(pv, lat, lon, pv_thresholds=None, n_thresholds=100)` | PV → φ_e table, ROI method. Same return format |
| `grid_cell_areas(lat, lon)` | Cell areas (m²) of a regular lat-lon grid, shape `(nlat, nlon)` |

If `pv_thresholds` is `None`, `n_thresholds` values are spaced linearly between the
minimum and maximum of the valid PV values.

### `eqlat.piecewise` / `eqlat.roi_fast`

| Function | Description |
|---|---|
| `eqlat_field_piecewise(pv, lat, lon, n_thresholds=200)` | φ_e at every grid point, shape `(nlat, nlon)` |
| `eqlat_field_roi(pv, lat, lon, n_thresholds=200)` | Same, with the ROI method |

### `eqlat.utils`

| Name | Description |
|---|---|
| `R_EARTH` | Earth radius, 6.371·10⁶ m |
| `total_sphere_area()` | 4πR² |
| `area_to_eqlat(area)` | A → φ_e = arcsin(1 − A / 2πR²) |
| `eqlat_to_area(eqlat)` | φ_e → A = 2πR² (1 − sin φ_e) |

### `eqlat.interpolation`

| Function | Description |
|---|---|
| `potential_temperature(T, p, p0=1000, kappa=0.286)` | θ = T (p₀/p)^κ, with p in hPa |
| `interpolate_to_theta(field, T, p_levels, theta_target, log_p=True)` | Column-wise interpolation of a `(nlev, nlat, nlon)` field to one θ surface. Also returns the pressure of that surface |
| `interpolate_to_theta_vectorized(field, T, p_levels, theta_target)` | Faster vectorised version, used by `batch` |

### `eqlat.batch`

| Function | Description |
|---|---|
| `process_netcdf(path, *, method, n_thresholds, n_workers, ...)` | φ_e for a file on θ levels |
| `process_pressure_netcdf(path, theta_levels, *, method, include_o3, p0, kappa, ...)` | Interpolation to θ + φ_e for a file on pressure levels |
| `process_{netcdf,pressure_netcdf}_{roi,piecewise}` | Shortcuts with `method` preset |

---

## Package layout

```
eqlat/
├── pyproject.toml
├── README.md
├── tests/
│   └── test_eqlat.py
└── src/eqlat/
    ├── __init__.py          public API, __version__
    ├── utils.py             grid cell areas, area ↔ φ_e conversion
    ├── piecewise.py         piecewise-constant method
    ├── roi_fast.py          ROI method (contourpy)
    ├── interpolation.py     θ computation, pressure → isentrope interpolation
    ├── batch.py             NetCDF batch processing (xarray)
    ├── calc/
    │   ├── process_fields.py    CLI: φ_e for a directory of reanalysis files
    │   └── merge_reanalysis.py  CLI: merge daily files
    └── download/
        ├── download_ERA_5.py
        ├── download_MERRA_2_new.py
        ├── download_mls.py
        └── download_ozonesondes.py
```

---

## Testing

```bash
pip install -e ".[dev]"
pytest
```

The tests check that the grid cell areas add up to the surface of the sphere, and that
both methods return the geographic latitude (within 1°) for a zonally symmetric PV
field.

---

## Known limitations

- **Regular grids only.** `grid_cell_areas` assumes evenly spaced latitudes and
  longitudes. Gaussian or reduced grids must be regridded first.
- **Non-monotonic θ profiles.** `interpolate_to_theta_vectorized` assumes θ is roughly
  monotonic in the column. Where a column crosses the target θ more than once (for
  example in a near-surface inversion), the crossing that comes last in the file's
  level order is used. The pressure levels are not sorted beforehand.
- **Hard-coded settings in `eqlat-process-fields`.** The θ levels and the file name
  patterns (`*pv*` / `*EPV*`, date taken from the file name) are fixed in the code.
- **MLS search limit.** `download_mls` asks `earthaccess` for at most 10 granules per
  product. Longer date ranges are truncated; split them into shorter ranges.
- **Ozonesonde `--url`** is effectively required, because there is no default station.

---

## License

Released under the [MIT License](LICENSE).

The ROI method in `roi_fast.py` is an independent implementation based on the algorithm
described in Añel et al. (2013). It is not derived from the authors' original IDL code
(GPL-3.0).
