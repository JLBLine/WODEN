
SRCLIST=srclist_singlepoint_power.yaml
FINE_CHAN_RES=40
NCHANS=32

METAFITS=../metafits/1126115208_metafits.fits

hyperdrive vis-simulate \
    -s $SRCLIST -m $METAFITS \
    --freq-res $FINE_CHAN_RES --num-fine-channels $NCHANS \
    --num-timesteps 14 \
    -o hyperdrive_dipflags.uvfits

hyperdrive vis-simulate \
    -s $SRCLIST -m $METAFITS \
    --freq-res $FINE_CHAN_RES --num-fine-channels $NCHANS \
    --num-timesteps 14 --unity-dipole-gains \
    -o hyperdrive_noflags.uvfits


METAFITS=../metafits/1088285600_DipAmps.metafits

hyperdrive vis-simulate \
    -s $SRCLIST -m $METAFITS \
    --freq-res $FINE_CHAN_RES --num-fine-channels $NCHANS \
    --num-timesteps 14 \
    -o hyperdrive_dipamps.uvfits

hyperdrive vis-simulate \
    -s $SRCLIST -m $METAFITS \
    --freq-res $FINE_CHAN_RES --num-fine-channels $NCHANS \
    --num-timesteps 14 --unity-dipole-gains \
    -o hyperdrive_perfect.uvfits