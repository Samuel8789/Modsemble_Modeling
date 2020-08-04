# !/bin/sh

echo "Do you want to install beta"
select yn in "Yes" "No"; do
	case $yn in
		Yes ) make; matlab -nodesktop -nosplash -batch "startup0;exit"; break;; 
		No ) exit
	esac
done

echo "Designated MATLAB MEX Directory"

make

echo "Installation Complete"




