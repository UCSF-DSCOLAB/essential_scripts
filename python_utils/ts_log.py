from datetime import datetime

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
This version does NOT actively flush the stdout buffer, so is meant for
interactive work or scripting that will be run with python's '-u' option
(e.g. python -u <script>) which eliminates the need for stdout flushing.
If you notice that your log messages are being held back until job
completion and cannot add this '-u' option in running the script, source
the 'ts_log__flush.py' version instead.
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
