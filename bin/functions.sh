#!/bin/bash

set -e
set -u
set -o pipefail

NBR='[+-]?[0-9]+\.?[0-9]*([eE][+-]?[0-9]+)?'
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
#
# DEPRECATED FUNCTION: USE $BIN/find_max instead!
#
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
# max()
# {
#     [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
#                           return 1; }
#     declare -r -i N="$1"
# 
# Sort according to -Im[eps^-1(omega)]
#     cat \
#         | "$BIN/convert" --from bin --to text --type "$TYPE" \
#         | nl \
#         | sed -E "s/$CMPL/\1\t\3/g" \
#         | awk -F'\t' '{ 
#               printf( "%i\t%.20E\t%.20E\t%.20E\n" \
#                     , $1, $2, $3, $3/($2**2 + $3**2) ); 
#           }' \
#         | sort -g -k4 \
#         | tail -n $N \
#         | tac
# }


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
#   freq1    n1    Re[eps_n1]    Im[eps_n1]    -Im[eps_n1^(-1)]    IPR[n1]
#   freq2    n2    Re[eps_n2]    Im[eps_n2]    -Im[eps_n2^(-1)]    IPR[n2]
#   freq3    n3    Re[eps_n3]    Im[eps_n3]    -Im[eps_n3^(-1)]    IPR[n3]
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
    echo -e "freq\tindex\treal\timag\tloss\tipr" >> "$spectrum"

    declare freq
    declare data
    declare loss
    declare ipr

    {
        echo -n "[*] MAX + IPR: " 1>&2
        tput sc 1>&2 # save cursor position
        for f in ${EPS_BASE}.*.eigenvalues.bin; do
            freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.eigenvalues.bin/\1/")

            tput rc 1>&2 # restore cursor position
            echo -n " frequency = $freq ..." 1>&2
            
            data=$( cat "$f" \
                  | $BIN/find_max --type $TYPE --number $N )
            loss=$( echo -e "$data" \
                  | awk '{ printf("%.20E", $3/($2**2 + $3**2)) }' )
            ipr=$( cat "${EPS_BASE}.${freq}.eigenstates.bin" \
                 | $BIN/ipr --type $TYPE --index "$(echo -e "$data" | cut -f 1)" )

            echo -e "$freq\t$data\t$loss\t$ipr"
        done 
        echo 1>&2
    } | sort -s -g -k 1 >> "$spectrum"
}


excitations()
{
    echo -n "[*] Calculating excitations: " 1>&2
    tput sc 1>&2 # save cursor position
    for f in ${EPS_BASE}.*.eigenvalues.bin; do
        freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.eigenvalues.bin/\1/")

        tput rc 1>&2 # restore cursor position
        echo -n " freq = $freq ..." 1>&2
        
        $BIN/reorder --type $TYPE \
                     --in.eigenvalues "${EPS_BASE}.${freq}.eigenvalues.bin" \
                     --in.eigenstates "${EPS_BASE}.${freq}.eigenstates.bin" \
                     --out.eigenvalues "reordered-${EPS_BASE}.${freq}.eigenvalues.bin" \
                     --out.eigenstates "reordered-${EPS_BASE}.${freq}.eigenstates.bin"
    done 
    echo 1>&2

}


reorder_all()
{
    declare freq

    echo -n "[*] Reordering: " 1>&2
    tput sc 1>&2 # save cursor position
    for f in ${EPS_BASE}.*.eigenvalues.bin; do
        freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.eigenvalues.bin/\1/")

        tput rc 1>&2 # restore cursor position
        echo -n " freq = $freq ..." 1>&2
        
        $BIN/reorder --type $TYPE \
                     --in.eigenvalues "${EPS_BASE}.${freq}.eigenvalues.bin" \
                     --in.eigenstates "${EPS_BASE}.${freq}.eigenstates.bin" \
                     --out.eigenvalues "reordered-${EPS_BASE}.${freq}.eigenvalues.bin" \
                     --out.eigenstates "reordered-${EPS_BASE}.${freq}.eigenstates.bin"
    done 
    echo 1>&2
}


loss_spectrum()
{
    declare -r direction=$( echo "$DIRECTION" \
                          | sed -E 's/\(\s*('$NBR')\s*,\s*('$NBR')\s*,\s*('$NBR')\s*\)/(\1,\3,\5)/g' )
    declare -r qs=$( \
python3 <<-EOF
begin = $Q_MIN
end   = $Q_MAX
step  = $Q_STEP
acc   = []
while begin <= end:
    acc.append("%.5E" % begin)
    begin += step
print(','.join(acc))
EOF
)
    echo "$qs" | tr ',' '\n' | wc -l >  ".temp.loss_spectrum.qs.dat"
    echo "$qs" | tr ',' '\n'         >> ".temp.loss_spectrum.qs.dat"

    declare freq
    declare spectrum

    tput sc 1>&2 # save cursor position
    for f in ${EPS_BASE}.*.matrix.bin; do
        freq=$(echo "$f" | sed -E "s/$EPS_BASE\.($NBR)\.matrix.bin/\1/")
        tput rc 1>&2 # restore cursor position
        echo -n "[*] frequency = $freq ..." 1>&2

        spectrum="${SPECTRUM/%.dat/.${freq}.${direction}.dat}"
        echo -e "# epsilon(q) for EPS_BASE = '$EPS_BASE'" > "$spectrum"
        echo -e "# freq = $freq" >> "$spectrum"
        echo -e "q\teps_r\teps_i" >> "$spectrum"
        $BIN/loss_function \
                --type "$TYPE" \
                --epsilon "$f" \
                --positions "$COORDINATES" \
                --direction "$direction" \
                --q "$qs" >> "$spectrum"
        cat  "$spectrum" \
            | tail -n +4 \
            | awk -F'\t' '{ printf("%.20E\n", $3 / ($2**2 + $3**2)); }' \
            | cat <(echo "$freq") - \
            > ".temp.loss_spectrum.spectrum.${freq}.dat"
    done
    echo 1>&2 # add a trailing newline

    echo "[*] Combining..." 1>&2
    paste .temp.loss_spectrum.qs.dat .temp.loss_spectrum.spectrum.*.dat \
        > "${SPECTRUM/%.dat/.${direction}.dat}"
    
    rm -f .temp.loss_spectrum.*.dat
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

    tput sc
    echo -n "[*] frequency = $freq, N = $N ..." 1>&2
    tput rc

    # Check whether we actually have to do anything
    [[ -f "$eigenmode" ]] && return 0

    echo -e "# Plasmon mode at $freq eV for EPS_BASE = '$EPS_BASE'" > "$eigenmode"
    echo -e "# max = $N" >> "$eigenmode"
    echo -e "x\ty\tz\treal\timag" >> "$eigenmode"

    declare -r -i index=$(cat "$spectrum" | grep -E "^${freq}" | cut -f 2)
    cat "${EPS_BASE}.${freq}.eigenstates.bin" \
        | "$BIN/convert" --type "$TYPE" --from bin --to text --column $index \
        | sed -E "s/$CMPL/\1\t\3/g" \
        | paste "$coordinates" "-" \
        | sort -s -g -k 1,2 \
        | awk ' NR == 1 { x = $1; } 
                NR > 1 && x != $1 { x = $1; printf("\n"); }
                { printf("%s\n", $0); } ' \
        >> "$eigenmode"
}


eigenstate_test()
{
    [[ "$#" -eq 2 ]] || { echo "$FUNCNAME takes exactly two arguments: " \
                               "<freq> and <N>." 1>&2
                          return -1; }
    # We want to store freq in the same "format" as the
    # 'plasmon_spectrum' function.
    declare -r freq=$(printf "%.6f" "$1")
    declare -r -i N="$2"
    declare -r coordinates="$COORDINATES"
    declare eigenmode

    tput sc 1>&2
    for i in $(seq $N); do
        tput rc 1>&2
        echo -n "[*] frequency = $freq, i = $i ..." 1>&2

        eigenmode="${EPS_BASE}.${freq}.${i}.eigenstate.dat"
        echo -e "# Plasmon mode at $freq eV for EPS_BASE = '$EPS_BASE'" > "$eigenmode"
        echo -e "# index = $i" >> "$eigenmode"
        echo -e "x\ty\tz\treal\timag" >> "$eigenmode"

        cat "${EPS_BASE}.${freq}.eigenstates.bin" \
            | "$BIN/convert" --type "$TYPE" --from bin --to text --column $i \
            | sed -E "s/$CMPL/\1\t\3/g" \
            | paste "$coordinates" "-" \
            | sort -s -g -k 1,2 \
            | awk ' NR == 1 { x = $1; } 
                    NR > 1 && x != $1 { x = $1; printf("\n"); }
                    { printf("%s\n", $0); } ' \
            >> "$eigenmode"
    done
    echo 1>&2
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


correlation_all()
{
    [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
                          return 1; }
    declare -r -i N="$1"
    declare -r spectrum="${SPECTRUM/%.dat/.${N}.dat}"
    declare -r corrfile="${CORRELATION/%.dat/.${N}.dat}"

    declare freq_1
    declare freq_2
    declare line_1
    declare line_2
    declare dot

    echo -e "# Correlation for N = $N" > "$corrfile"
    echo -e "freq\tcorr_r\tcorr_i" >> "$corrfile"

    {
        tput sc 1>&2
        IFS='' read -r line_1
        while IFS='' read -r line_2; do
            freq_1=$(echo -e "$line_1" | cut -f 1)
            freq_2=$(echo -e "$line_2" | cut -f 1)

            tput rc 1>&2
            echo -n "[*] freq_1 = ${freq_1}, freq_2 = ${freq_2} ..." 1>&2

            dot=$(correlation $freq_1 $freq_2 $N)

            echo -e "$freq_1\t$dot"
            line_1="$line_2"
        done 
        echo 1>&2
    } < <( cat "$spectrum" \
               | tail -n "+4" ) \
      >> "$corrfile"
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


ipr_all()
{
    [[ "$#" -eq 1 ]] || { echo "$FUNCNAME takes exactly one argument: <N>." 1>&2
                          return 1; }
    declare -r -i N="$1"
    declare -r spectrum="${SPECTRUM/%.dat/.${N}.dat}"
    declare -r iprfile="${IPR/%.dat/.${N}.dat}"

    declare freq
    declare line
    declare ipr

    echo -e "# IPR for N = $N" > "$iprfile"
    echo -e "freq\tipr" >> "$iprfile"

    tput sc 1>&2
    while IFS='' read -r line; do
        freq=$(echo -e "$line" | cut -f 1)

        tput rc 1>&2
        echo -n "[*] freq = ${freq} ..." 1>&2

        ipr=$(compute_ipr "$freq" "$N")

        echo -e "$freq\t$ipr"
    done < <( cat "$spectrum" \
               | tail -n "+4" ) \
         >> "$iprfile"
    echo 1>&2
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
                                  if(prev_sgn == -1 && cur_sgn != prev_sgn) {
                                      printf("%s\n%s\n", prev, cur);
                                  }
                                  prev = cur; prev_sgn = cur_sgn;
                              } ' ) \
         >> "$plasmons"
}



