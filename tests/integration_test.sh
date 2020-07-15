TROPOMI_DIR=/data/uptrop/nobackup/tropomi/Data/
PANDORA_DIR=/data/uptrop/nobackup/pandora/
OUT_DIR=./integration_outputs

echo Testing tropomi_ut_no2...
python ../uptrop/tropomi_ut_no2.py --season djf --grid_res 4x5 --cloud_product fresco --cloud_threshold 07 $TROPOMI_DIR $OUT_DIRi

echo Testing ut_no2_gc_test.py
python ../uptrop/ut_no2_gc_test.py --gc_dir $TROPOMI_DIR --out_path $OUT_DIR/gc_test.nc4

echo Testing compare_tropomi_pandora.py
python ../uptrop/compare_tropomi_pandora.py $TROPOMI_DIR $PANDORA_DIR $OUT_DIR --start_date 2019-06-01 --end_date 2019-06-30 --apply_bias_correction True

echo Testing fresco_cl_err.py
python ../uptrop/fresco_cld_err.py --s5p_data_dir $TROPOMI_DIR --output_dir $OUT_DIR
