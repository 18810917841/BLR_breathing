#!/usr/bin/expect -f
set user liyanrong
set host lxslc6.ihep.ac.cn
set password !qaz2wsx3edc
set file1 /home/yuyang/Desktop/Radiation_pressure/RAD_V6
set file2 /mbh/mbhd01/user/liyanrong/Yuyang/
set timeout -1
spawn scp -r $file1 $user@$host:$file2
expect "*assword:*"
send "$password\r"
expect eof
