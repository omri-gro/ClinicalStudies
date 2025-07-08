"""
SSH into given device (using a provided port number).
List all cases that were either scanned/analyzed or signed-off in the last N days.
Create directory in the device's /tmp directory.
For each case in list, get the event
"""
import paramiko


# move these to a utilities module later
def ssh_connect(hostname, username, password, port=22):
    # Create an SSH client
    client = paramiko.SSHClient()

    # Automatically add the server's host key
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        # Connect to the SSH server
        client.connect(hostname, port=port, username=username, password=password)

        return client
    except Exception as e:
        print(f"\033[91mError connecting to {hostname}: {e}\033[0m")
        return None


def download_file(sftp, remote_path, local_path):
    try:
        # Download the remote file
        sftp.get(remote_path, local_path)
        return True
    except Exception as e:
        print(f"\033[91mError downloading file {remote_path}: {e}\033[0m")
        return False


def is_remote_directory_exist(sftp, remote_path):
    try:
        # Attempt to change to the directory
        sftp.chdir(remote_path)
        return True
    except FileNotFoundError:
        # FileNotFoundError is raised if the directory doesn't exist
        return False
    except Exception as e:
        # Handle other exceptions (e.g., permission denied)
        print(f"Error checking directory existence: {e}")
        return False


def get_events_command(uuid, auth_key, host="localhost", dst="/tmp", save_name=None):
    """
    Create command for saving the events json of a case in a given location.
    Args:
        uuid (str): Case's  UUID
        auth_key (str): Such as swtoken
        host (str): Appears between 'https://' and '/analysis'
        dst (str): Directory to save the json in
        save_name (str): Name to save json as. If not provided will use <UUID>.json

    Returns:
        str: The full command to run
    """
    if not save_name:
        save_name = f'{uuid}.csv'

    command_str = f'curl GET "https://{host}/analysis/scans/{uuid}/_events" -H "accept: application/json" ' \
                  f'-H "Authorization: Bearer {auth_key}" --insecure --output {dst}/{save_name}.json'
    return command_str



# Remote device credentials
hostname = "remote_device_ip_or_hostname"
port = 22  # Default SSH port
username = "your_username"
password = "your_password"

# Commands to execute on the remote device
commands = [
    "ll -lrt",
    "curl -v --tlsv1.2 -X GET 'https://site_name/_events' -H 'accept: application/json' -H 'Authorization: Bearer long_key' --insecure"
]

try:
    # Create an SSH client
    ssh_client = paramiko.SSHClient()

    # Automatically add the remote server's SSH key (not recommended for production)
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # Connect to the remote server
    ssh_client.connect(hostname, port, username, password)
    print(f"Connected to {hostname}")

    # Execute commands
    for command in commands:
        print(f"Running command: {command}")
        stdin, stdout, stderr = ssh_client.exec_command(command)

        # Get command output and error
        output = stdout.read().decode()
        error = stderr.read().decode()

        # Print the output
        if output:
            print(f"Output:\n{output}")
        if error:
            print(f"Error:\n{error}")

    # Close the SSH connection
    ssh_client.close()
    print("Connection closed.")

except Exception as e:
    print(f"An error occurred: {e}")

