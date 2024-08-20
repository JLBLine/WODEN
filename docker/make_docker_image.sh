##Build some CUDA versions. Could we just do this with a multi-arch flag?
for arch in "60" "61" "70" "75" "80" "86"
do
    docker build --progress=plain --build-arg="CUDA_ARCH=${arch}" -t jlbline/woden-2.3:cuda-$arch -f Dockerfile_cuda .
    docker push jlbline/woden-2.3:cuda-$arch
done

##Build the setonix version. Someone could replace this with a HIP equivalent
##of the CUDA stuff above if they can work out how to do it.
docker build -t jlbline/woden-2.3:setonix -f Dockerfile_setonix .
docker push jlbline/woden-2.3:setonix


