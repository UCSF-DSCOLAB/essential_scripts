# Defines a simple function, ts_log, that prints a timestamped (ts) log message
# to stdout where the given 'msg' is printed after timestamp in the form of:
# [2025-03-11T14:05:47] hello
#
# Usage: ts_log hello
ts_log () {
    # Display help message if --help option is provided
    if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
        echo "Usage: ts_log [MESSAGE]"
        echo "prints a timestamped (ts) log message to stdout where the given"
        echo "'MESSAGE' is printed after a timestamp."
        echo ""
        echo "e.g. 'ts_log hello' will yield:"
        echo "[$(date +"%Y-%m-%dT%H:%M:%S")] hello"
    else
        echo "[$(date +"%Y-%m-%dT%H:%M:%S")] ""$*"
    fi
}
