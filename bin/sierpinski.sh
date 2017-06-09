#!/bin/bash
set -e
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
IPR=
CORRELATION=
COMMANDS=
DIRECTION=
Q_MIN=
Q_MAX=
Q_STEP=


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
            IPR="${i#*=}"
        ;;
        --correlation=*)
            CORRELATION="${i#*=}"
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
        --direction=*)
            DIRECTION="${i#*=}"
        ;;
        --q-min=*)
            Q_MIN="${i#*=}"
        ;;
        --q-max=*)
            Q_MAX="${i#*=}"
        ;;
        --q-step=*)
            Q_STEP="${i#*=}"
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


[[ "${TYPE}x" = "x" ]]        && TYPE="cdouble"
[[ "${SAMPLE_TYPE}x" = "x" ]] && SAMPLE_TYPE="sierpinski:carpet"
[[ "${SPECTRUM}x" = "x" ]]    && SPECTRUM="Spectrum.${START_WIDTH}.${DEPTH}.dat"
[[ "${PLASMONS}x" = "x" ]]    && PLASMONS="Plasmons.${START_WIDTH}.${DEPTH}.dat"
[[ "${IPR}x" = "x" ]]         && IPR="IPR.${START_WIDTH}.${DEPTH}.dat"
[[ "${CORRELATION}x" = "x" ]] && CORRELATION="Correlation.${START_WIDTH}.${DEPTH}.dat"
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
echo "[***] IPR              = ${IPR}" 1>&2
echo "[***] CORRELATION      = ${CORRELATION}" 1>&2
echo "[***] EPS_BASE         = ${EPS_BASE}" 1>&2
echo "[***] COORDINATES      = ${COORDINATES}" 1>&2
echo "[***] HAMILTONIAN      = ${HAMILTONIAN}" 1>&2
echo "[***] ENERGIES         = ${ENERGIES}" 1>&2
echo "[***] STATES           = ${STATES}" 1>&2
echo "[***] POTENTIAL        = ${POTENTIAL}" 1>&2
echo "[***] DIRECTION        = ${DIRECTION}" 1>&2
echo "[***] Q_MIN            = ${Q_MIN}" 1>&2
echo "[***] Q_MAX            = ${Q_MAX}" 1>&2
echo "[***] Q_STEP           = ${Q_STEP}" 1>&2
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
            [[ "${LATTICE_CONSTANT}x" = "x" ]] \
                && { echo "Need the '--lattice-constant' argument!" 1>&2; exit -1; }
            [[ "${HOPPING_VALUE}x" = "x" ]]    \
                && { echo "Need the '--hopping-value' argument!" 1>&2; exit -1; }
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
        ;;
        modes)
            # eigenstate_test 0.463500 512
            plasmon_modes
        ;;
        filter)
            filter_plasmons 1
        ;;
        loss)
            [[ "${DIRECTION}x" = "x" ]]  \
                && { echo "Need the '--direction' argument!" 1>&2; exit -1; }
            [[ "${Q_MIN}x" = "x" ]]      \
                && { echo "Need the '--q-min' argument!" 1>&2; exit -1; }
            [[ "${Q_MAX}x" = "x" ]]      \
                && { echo "Need the '--q-max' argument!" 1>&2; exit -1; }
            [[ "${Q_STEP}x" = "x" ]]     \
                && { echo "Need the '--q-tep' argument!" 1>&2; exit -1; }
            loss_spectrum
        ;;
        ipr)
            ipr_all 1 
        ;;
        correlation)
            correlation_all 1 
        ;;
        reorder)
            reorder_all
        ;;
        *)
            echo "Unknown command: $i" 1>&2
            exit -1
        ;;
    esac
done

