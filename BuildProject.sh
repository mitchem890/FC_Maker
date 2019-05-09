#Use this Docker Image to build a good base for you Docker recipt you can add things to this
#To add in more functionality
docker run --rm kaczmarj/neurodocker:master generate docker \
	--base ubuntu:16.04 --pkg-manager apt \
	--miniconda create_env=neuro \
		conda_install='python=3.6 scipy matplotlib numpy pandas' \
	--miniconda use_env=neuro \
		conda_install='jupyter' > Dockerfile



#Copy over my local files into the container
echo "COPY FCMaker /usr/" >> Dockerfile
echo "ENV PATH /opt/miniconda-latest/envs/neuro/bin:\$PATH" >> Dockerfile
#Set it to the 
echo "ENTRYPOINT [\"python\",\"-u\",\"/usr/FCDrawer.py\"]" >> Dockerfile

sed -i 's/apt-get/apt-get -y/g' Dockerfile
sed -i 's/nlibxmu-headers/libxmu-headers/g' Dockerfile
sed -i 's/nmesa-common-dev/mesa-common-dev/g' Dockerfile
sed -i '/libmng1/d' Dockerfile

docker build --tag fc_maker .
