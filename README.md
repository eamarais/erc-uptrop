# erc-uptrop Repository
TROPOMI cloud-slicing and validation code

A library for applying the cloud-slicing technique to TROPOMI total columns to obtain upper tropospheric (450-180 hPa) mixing ratios of NO2. This includes a feasibility test of the technique using synthetic observations from GEOS-Chem, validation and quantification of bias corrections by comparing TROPOMI and Pandora total columns at high-altitide sites, application of cloud-slicing using either the FRESCO-S or ROCINN-CAL cloud products, and comparison of effective cloud fractions and cloud top pressures from the two cloud products.

Examples of use
---------------

Cloud-slice synthetic partial columns from GEOS-Chem over North America at 4x5 in Spring:

```bash

   python ut_no2_gc_test.py --gc_dir "path/to/geoschem/data"  --out_path "/path/to/output/file" --resolution 4x5 --season mam --region NA
```
   
Co-sample daily Pandora and TROPOMI total NO2 columns at Izana in the first week of May:

```bash

   python compare_tropomi_pandora.py --trop_dir "/path/to/tropomi/data" --pan_dir "/path/to/pandora/data" --out_dir "/path/to/output/dir" --start_date 01-05-2019 --end_date 07-05-2019 --apply_bias_correction True --no2_col Trop --pandora_site izana
```
   
Cloud-slice Jun-Aug TROPOMI partial NO2 columns using the FRESCO-S cloud product in Summer:

```bash
   python tropomi_ut_no2.py --trop_dir "/path/to/tropomi/data/folder/" --out_dir "/path/to/output/folder/" --cloud_product fresco --season jja
```

Obtain January mean FRESCO-S and ROCINN-CAL cloud top pressures and fractions at 1x1:

```bash
   python fresco_cld_err.py  --trop_dir "/path/to/tropomi/data/folder/" --out_dir "/path/to/output/folder/" --start_date 01-01-2019 --end_date 31-01-2019 --out_res 1x1 --dlr_cld_top height
```

For a more advanced example bash script that runs each script over a given date range, see ``tests/integration_test.sh`` in your Uptrop folder.
See https://pages.github.io/erc-uptrop/docs/build/html/index.html?#scripts for detailed explanations of each function.

Documentation at https://pages.github.io/erc-uptrop/docs/build/html/index.html

Files copied from eam-group repo (https://github.com/eamarais/eam-group/tree/master/eloise/uptrop)

Citation where details of the method, and interpretation and visualization of the data is Marais et al., 2020.

ERC Starting Grant to address uncertainties in reactive nitrogen in the upper troposphere. 
Includes python code and routines to apply cloud-slicing to partial columns of NO2 to obtain mising ratios of NO2 in the upper tropposhepre. The routine
that does this is cloud_slice_ut_no2.py. 

