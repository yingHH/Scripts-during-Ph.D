# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-08-21 18:13:15
Last Modified by: Ying Huang
Last Modified time: 2020-08-21 18:13:15

Description: 
SSH port forwarding from local client to remote server.

Tip: Run `!Pyinstaller -F -w ./TransSSH.py` in IPython to pack this script to EXE file.
"""
import os
import sys
import time
import PySimpleGUI as sg
import wexpect


# global
#print = sg.Print

# funcs
def transSSH(password, user_name='yingh', ssh_IP='211.69.140.135', ssh_Port='2202', local_IP='localhost', local_Port='9000', remote_IP='211.69.140.135', remote_Port='9001', print=print):
    print('Run SSH')
    child = wexpect.spawn(
        'ssh -L {local_IP}:{local_Port}:{remote_IP}:{remote_Port} -p {ssh_Port} {user_name}@{ssh_IP}'.format(
            local_IP=local_IP,
            local_Port=local_Port,
            remote_IP=remote_IP,
            remote_Port=remote_Port,
            ssh_Port=ssh_Port,
            user_name=user_name,
            ssh_IP=ssh_IP,
        ),
        logfile = sys.stdout,
        encoding="utf-8",
    )

    print("""
    Opening SSH as follow information:

    user_name: {user_name}
    password: {password}

    ssh_IP: {ssh_IP}
    ssh_Port: {ssh_Port}

    remote_IP: {remote_IP}
    remote_Port: {remote_Port}

    local_IP: {local_IP}
    local_Port: {local_Port}
    """.format(
            password=password, user_name=user_name,
            ssh_IP=ssh_IP, ssh_Port=ssh_Port,
            remote_IP=remote_IP, remote_Port=remote_Port,
            local_IP=local_IP, local_Port=local_Port,
        ),
        background_color='skyblue',
    )
    #child.logfile = sys.stdout
    index = child.expect(['Permission Denied', 'password:'])

    sig = False

    if index == 1:
        try:
            child.sendline(password)
            print("Opened ssh", background_color='green')
            sig = True
        except:
            print("Wrong password:{}\n".format(password), background_color='red')
        #print(str(child))
    else:
        print("Error: \nPlease check your input as follow.\n{}".format(str(child)), background_color='red')
        
    return child, sig

def crt_wd1():
    layout = [
        [sg.Text("Information of server", font='Helvetica 16')],
        [sg.Text('IP')], [sg.InputText('211.69.140.135')], [sg.Text('Port')], [sg.InputText('2202')],
        [sg.Text('User Name')], [sg.InputText('yingh')],
        [sg.Text('Password')], [sg.InputText(password_char='*')],
        
        
        [sg.Text('Information of remote host', font='Helvetica 16')], 
        [sg.Text('IP')], [sg.InputText('211.69.140.135')], [sg.Text('Port')], [sg.InputText('9001')],
        
        [sg.Text('Information of local client', font='Helvetica 16')], 
        [sg.Text('IP')], [sg.InputText('localhost')], [sg.Text('Port')], [sg.InputText('9000')],
        
        [sg.Button('Submit')],
    ]
        
    return sg.Window(
        'Establish port forwarding through SSH from client',
        layout,
    )

def crt_wd2():
    layout = [
        [sg.MLine(key='-ML-'+sg.WRITE_ONLY_KEY, size=(70, 30))],
        [sg.Button('Back'), sg.Button('Refresh'), sg.Button('Close SSH process')],
    ]
    
    return sg.Window(
        'Runing outputs of SSH',
        layout,
        finalize=True,
    )

def crt_wd3():
    layout = [
        [sg.Text("Do you want to close SSH process?", key="-OUTPUT-")],
        [sg.Yes(), sg.No()],
    ]
    
    return sg.Window(
        'Close SSH checking',
        layout,
    )

# main

def main():
    # create GUI
    while True:
# Step 1
        wd1 = crt_wd1()
        
        event1, values1 = wd1.read()
        (
             ssh_IP, ssh_Port,
             user_name, password,
             remote_IP, remote_Port,
             local_IP, local_Port,
        ) = values1.values()

        if event1 == sg.WIN_CLOSED:
            wd1.close()
            break
        elif event1 == 'Submit':
            wd1.close()
# Step 2
        #is_transSSH_running = False

        while True:
            wd2 = crt_wd2()
            mprint = lambda *args, **kwargs: wd2['-ML-' + sg.WRITE_ONLY_KEY].print(*args, **kwargs)
            
            #if is_transSSH_running is False:
            process, sig = transSSH(
                password,
                user_name=user_name,
                ssh_IP=ssh_IP, ssh_Port=ssh_Port,
                remote_IP=remote_IP, remote_Port=remote_Port,
                local_IP=local_IP, local_Port=local_Port,
                print=mprint,            
            )

            if sig is True:
                mprint("Successfully established SSH port forwarding, holding windows ...", background_color='green')
            elif sig is False:
                mprint("Failed to establish SSH port forwarding, please check input information in blue", background_color='red')

            event2, values2 = wd2.read()

            if event2 == 'Back':
                process.close()
                wd2.close()
                break
            elif event2 == 'Refresh':
                process.close()
                wd2.close()
                continue
            elif event2 in (sg.WIN_CLOSED, 'Close SSH process'):
# Step 3
                wd3 = crt_wd3()
                event3, values3 = wd3.read()
                if event3 in ('Yes',):
                    wd3["-OUTPUT-"].update("Closing SSH ...")
                    process.close()
                    wd2.close()
                    wd3.close()
                    sg.Popup("SSH has been colsed.")
                    return # jump out the whole loop
                elif event3 in (sg.WIN_CLOSED, 'No'):
                    #is_transSSH_running = True
                    wd3.close()
                    continue

if __name__ == "__main__":
    #run_main()
    main()
    
