ZIP=transfer.zip
SRC=.
HOST=jxknight@niagara.scinet.utoronto.ca
DST=project/turnover
files=(
  @/code/config/
  @/code/epi-model/
  @/code/main/specs/
  @/code/main/system.py
  @/code/main/main.py
  @/code/main/variants.py
  @/code/main/batch.py
  @/code/main/compare.py
  @/code/main/surface.py
  @/code/main/surface-run.sh
  @/outputs/data/fit/
  @/outputs/data/surface/
)
fullfiles=(${files[@]/@/$SRC})     # expand pathnames
zip -r $ZIP ${fullfiles[@]}        # zip locally
scp $ZIP $HOST:$DST/$path          # copy zip file
ssh $HOST unzip $DST/$ZIP -d $DST  # unzip on scinet
ssh $HOST rm $DST/$ZIP             # remove the zipfile from scinet
rm $ZIP                            # remove the zipfile locally