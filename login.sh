#!/usr/bin/expect -f
set user liyanrong
set host lxslc6.ihep.ac.cn
set password !qaz2wsx3edc
set timeout 10
spawn ssh $user@$host
expect "*assword:*"
send "$password\r"
interact
