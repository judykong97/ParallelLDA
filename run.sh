# Read which interface to clean
process=1
iteration=1000
staleness=1
while test $# -gt 0; do
    case "$1" in 
        -p )
            process=$2
            shift 2
        ;;
        -i )
            iteration=$2
            shift 2
        ;;
        -s )
            staleness=$2
            shift 2
        ;;
        *)
            break
        ;;
    esac
done;

if [ $process == 1 ];
then
    ./lda $iteration $staleness
else 
    mpirun -np $process ./lda-mpi $iteration $staleness
fi
