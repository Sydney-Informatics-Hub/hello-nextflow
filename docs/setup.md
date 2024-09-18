# Setting up your computer

In this workshop, we will be using Virtual Machines (VM) on the
[ARDC Nectar Research Cloud](https://ardc.edu.au/services/ardc-nectar-research-cloud/).

The requirements for this workshop are a personal computer with:

- Visual Studio Code (VSCode)
- A web browser

Below, you will find instructions on how to set up VSCode and connect toyour VM.
Each participant will be provided with their instances IP address prior to the workshop.
Before the workshop, you must have the following:

1. VSCode installed.
2. The necessary VSCode extensions installed.
3. Be able to connect to your VM.  

## Installing Visual Studio Code

Visual Studio Code (VSCode) is a versatile code editor that we will use for the
workshop. We will use VSCode to connect to the VM, navigate the directories,
edit, view and download files. 

1. Download VSCode by following the 
[installation instructions](https://code.visualstudio.com/docs/setup/setup-overview)
for your local Operating System.  
2. Open VSCode to confirm it was installed correctly.  

![](img/vscode_0.png)

## Installing the VSCode extensions

Specific VSCode extensions are required to connect to the VM and make working
with Nextflow files easier (i.e. syntax highlighting).  

1. In the VSCode sidebar on the left, click on the extensions button (four
blocks).
2. In the Extensions Marketplace search bar, search for `remote ssh`. Select
**"Remote - SSH"**.

![](img/vscode_1.png)
3. Click on the blue `Install` button.

![](img/vscode_2.png)
4. Search for `nextflow` and install the **"Nextflow"** extension.  

![](img/vscode_3.png)
5. Close the Extensions tab and sidebar

## Setting up your remote SSH config  

1. In VSCode, press `Ctrl+Shift+P` (`Command+Shift+P` on mac) to open the Command Palette.  
 
![](img/ssh_0.png)
2. Type `remote ssh` and select **`Remote-SSH: Add New SSH Host...`**. This may
appear in a different position in the list.

![](img/ssh_1.png)
3. Enter the SSH connection string with the IP address that was provided to you. The connection string should look like **`ssh user1@XXX.XXX.XX.XX`**. Ensure that you replace the "XXX..." with your allocated IP address. Press `Enter`.

![](img/ssh_2.png)
4. You will be prompted to `Select SSH configuration file to update`. Select your `.ssh/config` file. 

![](img/ssh_3.png)
5. You should receive a pop-up informing that a host as been added!

## Connecting to the VM  

Ensure you have configured your SSH details.  

1. In VSCode, press `Ctrl+Shift+P` (`Command+Shift+P` on mac) to open the Command Palette.  
 
2. Type `remote ssh` and select **`Remote-SSH: Connect to Host...`**. This may
appear in a different position in the list.

![](img/vm_0.png)
3. Select the IP address that you have configured.  
4. A new VSCode window will open and prompt you for your password. Input your
allocated password and hit 'Enter'.

![](img/vm_1.png)
5. Once the blue square in the bottom-left of the VSCode window shows
 `SSH: XXX.XXX.XX.XX` - you have successfully connected to your instance!  

![](img/vm_2.png)  
6. Select the File Explorer on the left sidebar (icon with two pages) or press 
`Ctrl+Shift+E` (Mac: `Cmd+Shift+E`).  

![](img/vm_3.png)  
7. Select **`Open Folder`**

![](img/vm_4.png)  
8. The correct file path should be input by default (`/home/userX/`). Press 'OK'.
Finally, toggle the terminal in VSCode by pressing `Ctrl+j` (`Cmd+j` on mac).

> Add check for file explorer  
