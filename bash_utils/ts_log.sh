# Defines a simple function, ts_log, that prints a timestamped (ts) log message
# to stdout where the given 'msg' is printed after timestamp in the form of:
# [2025-03-11 14:05:47] hello
#
# Usage: ts_log hello
ts_log () {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] "$1
}
