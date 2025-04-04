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
'flush'             (Optional) boolean, controls whether to add a call to
                    'sys.stdout.flush()' after printing.  See note below.
*args, **kwargs     Strings, any item given to the function (not named as
                    'sep') will be concatenated together as the log message

Examples:
ts_log("hello world")
ts_log("hello", "world", sep = True, flush = True)

Additional details on 'flush' behavior:
Depending on how a script is run, python my hold its stdout in a buffer until
after the entire script finishes running. When 'flush = True', a call is added
to push the stdout buffer through early, ensuring your log-file recieves
these messages more quickly.
Alternatively, using python's `-u` option when running a script, e.g. 
'python -u <script>' alters the native behavior to where stdout is not held
until the script completes, eliminating the need for a stdout buffer flush.
    """
    msg = list()
    sep = ''
    do_flush = False
    if 'sep' in kwargs:
        sep = kwargs['sep']
        del kwargs['sep']
    if 'flush' in kwargs:
        do_flush = kwargs['flush']
        del kwargs['flush']
    if len(args)!=0:
        msg.extend(args)
    msg.extend(kwargs.values())
    print(f'[{datetime.now().strftime("%Y-%m-%dT%H:%M:%S")}] {sep.join([str(s) for s in msg])}')
    if do_flush:
        sys.stdout.flush()
