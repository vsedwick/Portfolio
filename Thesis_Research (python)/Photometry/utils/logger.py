import sys
import functools
import traceback
from io import StringIO
import builtins

def logger(logfile="log.txt"):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            log_output = StringIO()  # Capture printed output
            sys_stdout = sys.stdout  # Backup the original stdout
            sys_stderr = sys.stderr  # Backup the original stderr
            sys_input = input  # Backup the original input function

            try:
                sys.stdout = log_output  # Redirect stdout to log_output
                sys.stderr = log_output  # Redirect stderr to log_output
                
                # Redefine input to log user inputs
                def logged_input(prompt=""):
                    user_input = sys_input(prompt)
                    log_output.write(f"\nUser Input: {user_input}\n")
                    return user_input

                # Replace input with logged_input
                builtins.input = logged_input

                result = func(*args, **kwargs)
                return result

            except Exception as e:
                # Capture and log the error
                log_output.write("\nError Occurred:\n")
                log_output.write(traceback.format_exc())  # Log the full traceback
                raise e  # Optionally re-raise the exception

            finally:
                # Restore the original stdout, stderr, and input
                sys.stdout = sys_stdout
                sys.stderr = sys_stderr
                builtins.input = sys_input  # Restore original input

                # Write everything to the log file
                with open(logfile, "a") as f:
                    f.write(log_output.getvalue())
                    f.write("\n")  # Newline after each log

        return wrapper
    return decorator
