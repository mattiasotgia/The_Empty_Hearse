setup larsoft_data v1_02_02
source /exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/../localProducts_sbnana_v09_93_01_e26_prof/setup
mrbslp

refresh_token () {
    export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u)
    htgettoken -a htvaultprod.fnal.gov -i icarus
}

refresh_token
