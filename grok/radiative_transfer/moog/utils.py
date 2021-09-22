import os

import signal
import subprocess


# Get the path of MOOGSILENT/moogsilent.
for executable in ("MOOGSILENT", "moogsilent"):
    try:
        moogsilent_path = subprocess.check_output(
            f"which {executable}",
            shell=True
        )
    except subprocess.CalledProcessError:
        continue
    else:
        moogsilent_path = moogsilent_path.strip()
        acceptable_moog_return_codes = (0, )
#else:
#    raise RuntimeError(f"No moogsilent found on $PATH")

moogsilent_path = "/Users/arc/anaconda3/bin/MOOGSILENT"
acceptable_moog_return_codes = (0, )

    
def moogsilent(
        input_filename, 
        cwd=None, 
        timeout=30, 
        shell=False, 
        env=None,
        **kwargs
    ):
    """ 
    Execute a MOOGSILENT-compatible input file with a timeout after which it
    will be forcibly killed.

    :param input_filename:
        The full path of the MOOGSILENT-compatible input file.

    :param cwd: [optional]
        The current working directory to specify. If this is not specified, the
        directory of MOOGSILENT will be used.

    :param timeout: [optional]
        The number of seconds to wait before killing the process.

    :param shell: [optional]
        Whether to execute MOOGSILENT in a shell.

    :param env: [optional]
        A dictionary of environment variables to supply.
    """

    if cwd is None:
        cwd = os.path.dirname(input_filename)
    
    if env is None and len(os.path.dirname(moogsilent_path)) > 0:
        env = {"PATH": os.path.dirname(moogsilent_path)}

    p = subprocess.run(
        [moogsilent_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=cwd,
        input=os.path.basename(input_filename) + "\n"*100,
        encoding="ascii"
    )


    return (p.returncode, p.stdout, p.stderr)