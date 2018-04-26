# Read which interface to clean
process=1
while test $# -gt 0; do
    case "$1" in 
        -p )
            process=$2
            shift 2
        ;;
        *)
            break
        ;;
    esac
done;

if [ $process == 1 ];
then
    echo "sequential"
    ./lda
else 
    echo "parallel"
    mpirun -np $process ./lda
fi
