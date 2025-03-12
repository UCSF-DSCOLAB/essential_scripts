from datetime import datetime
import sys

def ts_log(*args, **kwargs):
    """
This function prints a timestamped (ts) log message where all elements
are printed after a timestamp in the form of:

[2025-03-11T14:05:47] {msg}

Arguments:
'sep'               (Optional) String, adjusts the separator used when
                    contenating multiple strings into the log message
*args, **kwargs     Strings, any item given to the function (not named as
                    'sep') will be concatenated together as the log message

Version note:
This version also actively flushes the stdout buffer, ensuring that you will
see the log line in your job log in a timely fashion.
When running a script with '-u' (e.g. python -u <script>), this buffer flushing
line is not needed, and the sourcing the 'ts_log.py' version is recommended.
    """
    msg = list()
    if 'sep' in kwargs:
        sep = kwargs['sep']
        del kwargs['sep']
    else:
        sep = ''
    if len(args)!=0:
        msg.extend(args)
    msg.extend(kwargs.values())
    print(f'[{datetime.now().strftime("%Y-%m-%dT%H:%M:%S")}] {sep.join([str(s) for s in msg])}')
    sys.stdout.flush()
