from datetime import datetime

def ts_log(msg):
    """
    This function prints a timestamped (ts) log message where the given 'msg'
    is printed after a timestamp in the form of:
    
    [2025-03-11 14:05:47] {msg}

    Arguments:
        'msg' - String, the desired text of the log message.
    
    Version note:
    This version does NOT actively flush the stdout buffer, so is meant for
    interactive work or scripting that will be run with python's '-u' option
    (e.g. python -u <script>) which eliminates the need for stdout flushing.
    If you notice that your log messages are being held back until job
    completion and cannot add this '-u' option in running the script, source
    the 'ts_log__flush.py' version instead.
    """
    print(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] {msg}')
