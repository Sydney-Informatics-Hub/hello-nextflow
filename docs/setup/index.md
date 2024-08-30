# Setup

## Before the workshop  

1. Install VSCode  
2. Install VSCode extensions  
3. Configure ssh key

### Step 1: Install Visual Studio Code (VSCode)

Visual Studio Code (VSCode) is a versatile code editor that we will use for the
workshop. Follow the 
[installation instructions](https://code.visualstudio.com/docs/setup/setup-overview)
for your local Operating System.  

### Step 2. Install VSCode extensions  

To work with Nextflow and connect to the virtual machine (VM), you will need to
install specific VSCode extensions.  

> Navigate to X; screenshot  

- `ms-vscode-remote.remote-ssh`
- `nextflow.nextflow`

### Step 3: Connect to the virtual machine (VM)  

> Screenshots  

In VSCode,  

1. Press `Ctrl+Shift+P` (`Command+Shift+P` on mac) to open the Command Palette.
2. Type `Remote-SSH: Connect to Host...` and select it.
3. Enter the SSH connection string that was provided to you. (`ssh <user>@<ip>`)

VSCode will open a new window connected to your VM.  

### Step 4: Validate your setup  

> screenshots  

> Run a test Nextflow workflow  

> Run docker without root access
