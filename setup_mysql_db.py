from ebrowse import PATHS
import getpass
import os

### The commands are taken from: https://dzone.com/articles/setting-up-mysql-without-root-access

def mk_paths(path):
    if not os.path.exists(path):
        os.makedirs(path)

base_dir = PATHS['mysqldb_basepath']
data_dir = base_dir + "/data"
log_dir = base_dir + "/log"
tmp_dir = base_dir + "/tmp"
user_name = getpass.getuser()
port = PATHS['mysqldb_port']
mysql_share_dir = PATHS['mysql_share_dir']
mysql_scripts_dir = PATHS['mysql_scripts_dir']
mysql_install_dir = PATHS['mysql_base_dir']
socket = base_dir + "/socket"

mysql_cnf = """
# For explanations see
# http://dev.mysql.com/doc/mysql/en/server-system-variables.html

# This will be passed to all mysql clients
# It has been reported that passwords should be enclosed with ticks/quotes
# escpecially if they contain "#" chars...
# Remember to edit /etc/mysql/debian.cnf when changing the socket location.
[client]
port = {port}
socket = {base_dir}/socket

# Here is entries for some specific programs
# The following values assume you have at least 32M ram

# This was formally known as [safe_mysqld]. Both versions are currently parsed.
[mysqld_safe]
socket		= {base_dir}/socket
nice		= 0
log-error   = {log_path}/errlog.log
general_log_file   = {log_path}/genlog.log
pid-file	= {base_dir}/mysql.pid

[mysqld]
#
# * Basic Settings
#
user		= {user}
pid-file	= {base_dir}/mysql.pid
socket		= {base_dir}/socket
port		= {port}
basedir		= {mysql_install_dir}
datadir		= {data_dir}
tmpdir		= {tmp_dir}
log-error   = {log_path}/errlog.log
general_log_file   = {log_path}/genlog.log
#lc-messages-dir= {mysql_share_dir}
skip-external-locking
skip-character-set-client-handshake
default-storage-engine = InnoDB
character-set-server = utf8
transaction-isolation = READ-COMMITTED
""".format(user= user_name, mysql_share_dir=mysql_share_dir, base_dir = base_dir,
           port = str(port), data_dir=data_dir, tmp_dir = tmp_dir, mysql_install_dir= mysql_install_dir,
           log_path = log_dir)

initialise_mysql = "mysqld --user={user} --datadir={data_dir} --basedir={base_dir} --log-error={log_dir}/mysql.err " \
                   "--pid-file={base_dir}/mysql.pid  --socket={base_dir}/socket --port={port} --lc-messages-dir {base_dir}/share " \
                   "--initialize".format(user=user_name, base_dir = base_dir, data_dir =data_dir, log_dir = log_dir, port = str(port))

initialise_mysql = "mysqld --defaults-extra-file={conf}".format(conf=PATHS['mysqldb_basepath'] + "/my.cnf")

install_db_mysql = "{mysql_scripts_dir}/mysql_install_db --basedir={basedir} --datadir={datadir}".format(
    mysql_scripts_dir=mysql_scripts_dir, basedir = base_dir, datadir=data_dir)


mk_paths(base_dir)
mk_paths(data_dir)
mk_paths(log_dir)
mk_paths(tmp_dir)
mk_paths(mysql_share_dir + "/share")

with open(PATHS['mysqldb_basepath'] + "/my.cnf", "w") as ofh:
    ofh.write(mysql_cnf)

#os.system(initialise_mysql)
#print(install_db_mysql)
#os.system(install_db_mysql)


# what actually worked:
#/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/scripts/mysql_install_db --basedir=/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser --datadir=/nfs/research1/stegle/users/rkreuzhu/webapp_data/mysql/data


"""
/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysqladmin -u root password 'ebrowseWHO1000'
/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysqladmin -u root -h ebi-cli-001.ebi.ac.uk password 'new-password'

Alternatively you can run:
/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysql_secure_installation

which will also give you the option of removing the test
databases and anonymous user created by default.  This is
strongly recommended for production servers.

See the manual for more instructions.

You can start the MySQL daemon with:
cd /nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser ; /nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysqld_safe &

You can test the MySQL daemon with mysql-test-run.pl
cd /nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/mysql-test ; perl mysql-test-run.pl
"""

#/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysqld_safe --defaults-extra-file=/nfs/research1/stegle/users/rkreuzhu/webapp_data/mysql/my.cnf &
#/nfs/research1/stegle/users/rkreuzhu/conda-envs/ebrowser/bin/mysqladmin --port 3306 --socket /nfs/research1/stegle/users/rkreuzhu/webapp_data/mysql/socket -u root password 'ebrowseWHO100'

# then setup another user:
#mysql --socket /nfs/research1/stegle/users/rkreuzhu/webapp_data/mysql/socket --port 3306 -u root -p
# then do:
"GRANT ALL PRIVILEGES ON *.* TO 'ebrowse'@'localhost' IDENTIFIED BY 'mergedEqtlData';"
# create the database:
"CREATE DATABASE eqtldata;"
# and:
"quit"




# to shut the server down do
# mysqladmin --socket=/nfs/research1/stegle/users/rkreuzhu/webapp_data/mysql/socket shutdown -u root -p

"""
export BASE=/srv/tempo/pgms/mysql
$BASE/bin/mysqladmin 
--socket=$BASE/socket shutdown -u root -p
"""


