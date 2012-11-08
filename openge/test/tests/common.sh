cd $(dirname $0)

OGE="../../../../bin/openge"
DATA="../../data"

echo "OGE: $OGE"
echo `ls $OGE`
function err() {
        echo
        echo "ERROR: $1"
        echo
        exit 1
}
