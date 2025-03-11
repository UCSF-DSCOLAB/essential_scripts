from datetime import datetime
import sys

def ts_log(msg):
    """
    This function prints a timestamped (ts) log message where the given 'msg'
    is printed after a timestamp in the form of:
    
    [2025-03-11 14:05:47] {msg}

    Arguments:
        'msg' - String, the desired text of the log message.
    
    Version note:
    This version also actively flushes the stdout buffer, ensuring that you
    will see the log line in your job log in a timely fashion.
    When running a script with '-u' (e.g. python -u <script>), this buffer
    flushing line is not needed, and the sourcing the 'ts_log.py' version
    is recommended.
    """
    print(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] {msg}')
    sys.stdout.flush()
