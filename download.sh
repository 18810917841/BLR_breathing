#!/usr/bin/expect -f
set user liyanrong
set host lxslc6.ihep.ac.cn
set password !qaz2wsx3edc
set file1 /mbh/mbhd01/user/liyanrong/Yuyang/RAD_V6/data/*
set file2 /home/yuyang/Desktop/Radiation_pressure/RAD_V6/data
set timeout -1
spawn scp -r $user@$host:$file1 $file2
expect "*assword:*"
send "$password\r"
expect eof
