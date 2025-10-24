
woden_version="2.7"

#Build some CUDA versions.

for arch in "60" "61" "70" "75" "80" "86" "60;61;70;75;80;86"
do
    if [[ "$arch" == "60;61;70;75;80;86" ]]; then
        tag=multi
    else
        tag=${arch}
    fi

    docker build --no-cache --progress=plain --build-arg="CUDA_ARCH=${arch}" --build-arg="USE_BUILD=production" \
                 -t jlbline/woden-${woden_version}:cuda-${tag} -f Dockerfile_cuda .
    docker push jlbline/woden-${woden_version}:cuda-${tag}

    # # Uncomment the following lines if you want to build debug versions as well.
    # docker build --no-cache --progress=plain --build-arg="CUDA_ARCH=${arch}" --build-arg="USE_BUILD=debug" \
    #              -t jlbline/woden-${woden_version}:cuda-${tag}_debug -f Dockerfile_cuda .
    # docker push jlbline/woden-${woden_version}:cuda-${tag}_debug
done

##Build the setonix version. Someone could replace this with a HIP equivalent
##of the CUDA stuff above if they can work out how to do it.
docker build --no-cache --progress=plain --build-arg="USE_BUILD=production" \
    -t jlbline/woden-${woden_version}:setonix -f Dockerfile_setonix .
docker push jlbline/woden-${woden_version}:setonix

##If you need to experiment with environments, via a sandbox, you can do it like this:
# apptainer build --sandbox sandbox_ubuntu docker://nvidia/cuda:12.5.1-devel-ubuntu24.04

##Then when running. you need to run sudo outside the box to be sudo inside the box:
# sudo apptainer shell --writable sandbox_ubuntu

