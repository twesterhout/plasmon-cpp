#!/bin/bash

set -e
set -u
set -o pipefail

NBR='[+-]?[0-9]+\.?[0-9]*([eE][+-][0-9]+)?'
CMPL='\(('$NBR'),('$NBR')\)'


###############################################################################
# Generates the system if it doesn't exist.
# Calls the do_generate to perform all the hard work
#
# Variables:
# * [global] HAMILTONIAN   -- file name where tipsi stores H
# * [global] COORDINATES   -- file name where tipsi stores (X,Y,Z)
#
# Functions:
# * do_generate -- does the actual work. It should also use HAMILTONIAN and
#                  COORDINATES variables as filenames.
###############################################################################
generate()
{
    declare -r H_text="${HAMILTONIAN/%.bin/%.dat}"
    declare -r X_text="${COORDINATES/%.bin/%.dat}"

    if [[ -f "$H_text" && -f "$X_text" ]]; then
        # We check that the number of lines containing the actual data,
        # that is not starting with a '#', is the same for both files.
        [[    $(cat "$H_text" | grep -E -v '^\s*#' | wc -l) \
          -eq $(cat "$X_text" | grep -E -v '^\s*#' | wc -l) \
        ]] || { echo "[-] Hamiltonian and Coordinates files contain" \
                     "different number of points." 1>&2
                return -1; }
        echo "[+] Systems already exists." 1>&2
    else
        echo "[*] Generating the system..." 1>&2
        do_generate
        echo "[+] Successfully generated the system." 1>&2
    fi 
}


###############################################################################
# Converts the Hamiltonian from its text version (produced by the 
# 'generate_system.py' script) to binary archive readable with boost::archive.
#
# Variables
# * BIN           -- location of binaries
# * HAMILTONIAN   -- file name where tipsi stores H
# * TYPE          -- type of elements in H
###############################################################################
convert_hamiltonian()
{
    declare -r H_text="${HAMILTONIAN/%.bin/.dat}"
    declare -r H_bin="${HAMILTONIAN/%.dat/.bin}"

    if [[ -f "$H_bin" ]]; then
        echo "[+] Binary version of the Hamiltonian already exists." 1>&2
    else
        [[ -f "$H_text" ]] || { echo "[-] Hamiltonian not found."; return -1; }
        echo "[*] Converting: $H_text --> $H_bin ..." 1>&2
        cat "$H_text" | grep -E -v '^\s*#' \
                      | "$BIN/convert" --type "$TYPE" \
                                       --from "text" \
                                       --to "bin" \
                      > "$H_bin"
        echo "[+] Successfully converted the Hamiltonian." 1>&2
    fi
}


###############################################################################
# Calculates the Potential from atomic sites coordinates.
#
# Variables:
# * BIN           -- location of binaries
# * POTENTIAL     -- file name where to store V (in binary form).
# * COORDINATES   -- file name where tipsi stores {(x,y,z)}_i
# * TYPE          -- type of elements in V
###############################################################################
calculate_potential()
{
    declare -r V_bin="${POTENTIAL}"
    declare -r X_text="${COORDINATES}"

    if [[ -f "$V_bin" ]]; then
        echo "[+] Potential already exists." 1>&2
    else
        echo "[*] Calculating Potential ..." 1>&2
        cat "$X_text" | grep -E -v '^\s*#' \
                      | "$BIN/potential" --type "$TYPE" > "$V_bin"
        echo "[+] Successfully calculated the Potential." 1>&2
    fi
}


###############################################################################
# Given all eigenvalues of the dielectric function, prints some useful
# information about the first <N> maxima, where the values are sorted by
# -Im[eps_n^(-1)]. Eigenvalues are read as a binary archive from stdin.
# 
# This function takes one argument: <N>, i.e. the number of maxima to 
# print.
#
# Variables:
# * CMPL          -- regex for a complex number (defined inside this file).
# * BIN           -- location of the binaries.
# * TYPE          -- type of eigenvalues of Epsilon.
#
# Result:
#   n1    Re[eps_n1]    Im[eps_n1]    -Im[eps_n1^(-1)]
#   n2    Re[eps_n2]    Im[eps_n2]    -Im[eps_n2^(-1)]
#   ...
#   nN    Re[eps_nN]    Im[eps_nN]    -Im[eps_nN^(-1)]
###############################################################################
max()
{
    [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
                          return 1; }
    declare -r -i N="$1"

    cat \
        | "$BIN/convert" --from bin --to text --type "$TYPE" \
        | nl \
        | sed -E "s/$CMPL/\1\t\3/g" \
        | awk -F'\t' '{ 
              printf( "%i\t%.20E\t%.20E\t%.20E\n" \
                    , $1, $2, $3, $3/($2**2 + $3**2) ); 
          }' \
        | sort -g -k4 \
        | tail -n $N \
        | tac
}


###############################################################################
# Creates the "complete" plasmon spectrum, i.e. from all available 
# frequencies.
#
# Variables:
# * NBR           -- regex for a real number (defined in this file).
# * EPS_BASE      -- type of eigenvalues of Epsilon.
# * SPECTRUM      -- file name where to save the spectrum
# 
# Result:
#   freq1    n1    Re[eps_n1]    Im[eps_n1]    -Im[eps_n1^(-1)]
#   freq2    n2    Re[eps_n2]    Im[eps_n2]    -Im[eps_n2^(-1)]
#   freq3    n3    Re[eps_n3]    Im[eps_n3]    -Im[eps_n3^(-1)]
#   ...
#
###############################################################################
plasmon_spectrum()
{
    [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
                          return 1; }
    declare -r -i N="$1"
    declare -r spectrum="${SPECTRUM/%.dat/.${N}.dat}"

    echo -e "# Plasmon spectrum for EPS_BASE = '$EPS_BASE'" > "$spectrum"
    echo -e "# max = $N" >> "$spectrum"
    echo -e "freq\tindex\treal\timag\tloss" >> "$spectrum"

    declare freq
    for f in ${EPS_BASE}.*.eigenvalues.bin; do
        freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.eigenvalues.bin/\1/")
        echo "[*] frequency = $freq ..." 1>&2
        cat "$f" \
            | max $N \
            | tail -n 1 \
            | awk -F'\t' '{ printf("%.6f\t%s\n", '$freq', $0); }'
    done | sort -s -g -k 1 >> "$spectrum"
}


###############################################################################
# Extracts plasmon mode at given frequency
###############################################################################
eigenstate()
{
    [[ "$#" -eq 2 ]] || { echo "$FUNCNAME takes exactly two arguments: " \
                               "<freq> and <N>." 1>&2
                          return -1; }
    # We want to store freq in the same "format" as the
    # 'plasmon_spectrum' function.
    declare -r freq=$(printf "%.6f" "$1")
    declare -r -i N="$2"
    declare -r spectrum="${SPECTRUM/%.dat/.${N}.dat}"
    declare -r coordinates="$COORDINATES"
    declare -r eigenmode="${EPS_BASE}.${freq}.${N}.eigenstate.dat"

    echo "[*] frequency = $freq, N = $N ..." 1>&2

    # Check whether we actually have to do anything
    [[ -f "$eigenmode" ]] && return 0

    echo -e "# Plasmon mode at $freq eV for EPS_BASE = '$EPS_BASE'" > "$eigenmode"
    echo -e "# max = $N" >> "$eigenmode"
    echo -e "x\ty\tz\treal\timag" >> "$eigenmode"

    declare -r -i index=$(cat "$spectrum" | grep -E "^${freq}" | cut -f 2)
    cat "${EPS_BASE}.${freq}.eigenstates.bin" \
        | "$BIN/convert" --type "$TYPE" --from bin --to text \
        | cut -f $index \
        | sed -E "s/$CMPL/\1\t\3/g" \
        | paste "$coordinates" "-" \
        | sort -s -g -k 1,2 \
        | awk ' NR == 1 { x = $1; } 
                NR > 1 && x != $1 { x = $1; printf("\n"); }
                { printf("%s\n", $0); } ' \
        >> "$eigenmode"
}


plasmon_modes()
{
    declare freq

    for f in ${EPS_BASE}.*.eigenvalues.bin; do
        freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.eigenvalues.bin/\1/")
        eigenstate "$freq" 1
        eigenstate "$freq" 2
    done
}


correlation()
{
    [[ "$#" -eq 3 ]] || { echo "$FUNCNAME takes exactly three arguments:" \
                               "<freq_1>, <freq_2> and <N>." 1>&2
                          return -1; }
    declare -r freq_1=$(printf "%.6f" "$1")
    declare -r freq_2=$(printf "%.6f" "$2")
    declare -r -i N="$3"

    [[ $N -gt 0 && $N -lt 3 ]] || { echo "<max> may be either 1 or 2" 1>&2
                                            return 1; }

    paste <(cat "${EPS_BASE}.${freq_1}.${N}.eigenstate.dat" | cut -f 4,5) \
          <(cat "${EPS_BASE}.${freq_2}.${N}.eigenstate.dat" | cut -f 4,5) \
        | awk -F'\t' '
              BEGIN { acc[0] = 0; acc[1] = 0; } 
              { x[0] = $1; x[1] = $2;
                y[0] = $3; y[1] = $4;
                acc[0] += x[0] * y[0] + x[1] * y[1];
                acc[1] += x[0] * y[1] - x[1] * y[0];
              }
              END { printf("%.20E\t%.20E\n", acc[0], acc[1]); } '
}


compute_ipr()
{
    [[ "$#" -eq 2 ]] || { echo "$FUNCNAME takes exactly two arguments: " \
                               "<freq> and <N>." 1>&2
                          return 1; }
    declare -r freq=$(printf "%.6f" "$1")
    declare -r -i N="$2"

    cat "${EPS_BASE}.${freq}.${N}.eigenstate.dat" \
        | awk -F'\t' ' BEGIN { acc = 0; }
                       { ipr += ($4**2 + $5**2)**2; } 
                       END { printf("%.20E\n", ipr); } '
}


filter_plasmons()
{
    [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
                          return 1; }
    declare -r -i N="$1"
    declare -r spectrum="${SPECTRUM/%.dat/.${N}.dat}"
    declare -r plasmons="${PLASMONS/%.dat/.${N}.dat}"

    declare freq_1
    declare freq_2
    declare line_1
    declare line_2
    declare ipr_1
    declare ipr_2
    declare dot

    echo -e "# Plasmons for N = $N" > "$plasmons"
    echo -e "freq\tindex\treal\timag\tloss\tipr\tcorr_r\tcorr_i" >> "$plasmons"

    while IFS='' read -r line_1; do
        IFS='' read -r line_2

        freq_1=$(echo -e "$line_1" | cut -f 1)
        freq_2=$(echo -e "$line_2" | cut -f 1)

        echo "[*] freq_1 = ${freq_1}, freq_2 = ${freq_2} ..." 1>&2

        ipr_1=$(compute_ipr "$freq_1" "$N")
        ipr_2=$(compute_ipr "$freq_2" "$N")
        dot=$(correlation $freq_1 $freq_2 $N)

        echo -e "$line_1\t$ipr_1\t$dot"
        echo -e "$line_2\t$ipr_2\t$dot"
        echo
    done < <( cat "$spectrum" \
               | tail -n "+4" \
               | awk -F'\t' ' function sign(x) { return x >= 0 ? 1 : -1; }
                              NR == 1 { prev = $0; prev_sgn = sign($3); }
                              NR > 1 {
                                  cur = $0; cur_sgn = sign($3);
                                  if(cur_sgn != prev_sgn) {
                                      printf("%s\n%s\n", prev, cur);
                                  }
                                  prev = cur; prev_sgn = cur_sgn;
                              } ' ) \
         >> "$plasmons"
}



