#!/bin/bash
# set -e
set -u
# set -v
set -o pipefail

TYPE=
START_WIDTH=
DEPTH=
LATTICE_CONSTANT=
HOPPING_VALUE=
EPS_BASE=
COORDINATES=
HAMILTONIAN=
ENERGIES=
STATES=
POTENTIAL=
SAMPLE_TYPE=
SPECTRUM=
PLASMONS=
IPR_CORR=
COMMANDS=
WAVEVECTOR=


for i in "$@"; do
    case $i in
        --type=*)
            TYPE="${i#*=}"
        ;;
        --start-width=*)
            START_WIDTH="${i#*=}"
        ;;
        --depth=*)
            DEPTH="${i#*=}"
        ;;
        --sample-type=*)
            SAMPLE_TYPE="${i#*=}"
        ;;
        --lattice-constant=*)
            LATTICE_CONSTANT="${i#*=}"
        ;;
        --hopping-value=*)
            HOPPING_VALUE="${i#*=}"
        ;;
        --eps-base=*)
            EPS_BASE="${i#*=}"
        ;;
        --spectrum=*)
            SPECTRUM="${i#*=}"
        ;;
        --plasmon=*)
            PLASMONS="${i#*=}"
        ;;
        --ipr=*)
            IPR_CORR="${i#*=}"
        ;;
        --coordinates=*)
            COORDINATES="${i#*=}"
        ;;
        --hamiltonian=*)
            HAMILTONIAN="${i#*=}"
        ;;
        --energies=*)
            ENERGIES="${i#*=}"
        ;;
        --states=*)
            STATES="${i#*=}"
        ;;
        --potential=*)
            POTENTIAL="${i#*=}"
        ;;
        --wave-vector=*)
            WAVEVECTOR="${i#*=}"
        ;;
        --commands=*)
            COMMANDS="${i#*=}"
        ;;
        *)
            echo "[---] Unknown option: $i"
            exit -1
        ;;
    esac
done

[[ "${START_WIDTH}x" = "x" ]]      && { echo "Need the '--start-width' argument!" 1>&2; exit -1; }
[[ "${DEPTH}x" = "x" ]]            && { echo "Need the '--depth' argument!" 1>&2; exit -1; }
[[ "${LATTICE_CONSTANT}x" = "x" ]] && { echo "Need the '--lattice-constant' argument!" 1>&2; exit -1; }
[[ "${HOPPING_VALUE}x" = "x" ]]    && { echo "Need the '--hopping-value' argument!" 1>&2; exit -1; }
[[ "${WAVEVECTOR}x" = "x" ]]       && { echo "Need the '--wave-vector' argument!" 1>&2; exit -1; }
[[ "${COMMANDS}x" = "x" ]]         && { echo "Need the '--commands' argument!" 1>&2; exit -1; }
[[ "${TYPE}x" = "x" ]]        && TYPE="cdouble"
[[ "${SAMPLE_TYPE}x" = "x" ]] && SAMPLE_TYPE="sierpinski:carpet"
[[ "${SPECTRUM}x" = "x" ]]    && SPECTRUM="Spectrum.${START_WIDTH}.${DEPTH}.dat"
[[ "${PLASMONS}x" = "x" ]]    && PLASMONS="Plasmons.${START_WIDTH}.${DEPTH}.dat"
[[ "${IPR_CORR}x" = "x" ]]    && IPR_CORR="IPR.${START_WIDTH}.${DEPTH}.dat"
[[ "${EPS_BASE}x" = "x" ]]    && EPS_BASE="Epsilon.${START_WIDTH}.${DEPTH}"
[[ "${HAMILTONIAN}x" = "x" ]] && HAMILTONIAN="Hamiltonian.${START_WIDTH}.${DEPTH}.dat"
[[ "${COORDINATES}x" = "x" ]] && COORDINATES="Coordinates.${START_WIDTH}.${DEPTH}.dat"
[[ "${ENERGIES}x" = "x" ]]    && ENERGIES="Energies.${START_WIDTH}.${DEPTH}.bin"
[[ "${STATES}x" = "x" ]]      && STATES="States.${START_WIDTH}.${DEPTH}.bin"
[[ "${POTENTIAL}x" = "x" ]]   && POTENTIAL="Potential.${START_WIDTH}.${DEPTH}.bin"


echo "[***] TYPE             = ${TYPE}" 1>&2
echo "[***] START_WIDTH      = ${START_WIDTH}" 1>&2
echo "[***] DEPTH            = ${DEPTH}" 1>&2
echo "[***] LATTICE_CONSTANT = ${LATTICE_CONSTANT}" 1>&2
echo "[***] HOPPING_VALUE    = ${HOPPING_VALUE}" 1>&2
echo "[***] SAMPLE_TYPE      = ${SAMPLE_TYPE}" 1>&2
echo "[***] SPECTRUM         = ${SPECTRUM}" 1>&2
echo "[***] PLASMONS         = ${PLASMONS}" 1>&2
echo "[***] IPR_CORR         = ${IPR_CORR}" 1>&2
echo "[***] EPS_BASE         = ${EPS_BASE}" 1>&2
echo "[***] COORDINATES      = ${COORDINATES}" 1>&2
echo "[***] HAMILTONIAN      = ${HAMILTONIAN}" 1>&2
echo "[***] ENERGIES         = ${ENERGIES}" 1>&2
echo "[***] STATES           = ${STATES}" 1>&2
echo "[***] POTENTIAL        = ${POTENTIAL}" 1>&2
echo "[***] WAVEVECTOR       = ${WAVEVECTOR}" 1>&2
echo "[***] COMMANDS         = ${COMMANDS}" 1>&2


NBR='[+-]?[0-9]+\.?[0-9]*([eE][+-][0-9]+)?'
CMPL='\(('$NBR'),('$NBR')\)'


do_generate()
{
    ${BIN}/generate_system.py --type "${SAMPLE_TYPE}" \
                              --lattice-constant "${LATTICE_CONSTANT}" \
                              --hopping-value "${HOPPING_VALUE}" \
                              --start "${START_WIDTH}" \
                              --depth "${DEPTH}"
}

source $BIN/functions.sh

IFS=',' read -ra COMMANDS_LIST <<< "${COMMANDS}"
for i in "${COMMANDS_LIST[@]}"; do
    case $i in
        generate)
            generate
        ;;
        convert)
            convert_hamiltonian
        ;;
        potential)
            calculate_potential
        ;;
        spectrum)
            plasmon_spectrum 1
            plasmon_spectrum 2
        ;;
        modes)
            plasmon_modes
        ;;
        filter)
            filter_plasmons 1
        ;;
        loss)
            loss_spectrum "$WAVEVECTOR"
        ;;
        *)
            echo "Unknown command: $i" 1>&2
            exit -1
        ;;
    esac
done

