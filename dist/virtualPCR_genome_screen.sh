#!/usr/bin/env bash

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR_PATH="${SCRIPT_DIR}/virtualPCR.jar"
WORK_DIR="$(pwd)"

temp_input_dir=""
temp_java_output_dir=""
temp_run_dir=""
run_config=""
temp_java_output_path=""

show_help() {
	echo "Usage: $0 -i <genome1> [genome2 ...] -p <primers.fasta> [options]"
	echo ""	
	echo "Wrapper script for virtualPCR.jar. This script handles multiple genome inputs and generates a summary output file."
	echo "Make sure virtualPCR.jar is in the same directory as this script."
	echo ""
	echo "Required:"
	echo "  -i   One or more genome input files (supports multiple values and repeated -i)"
	echo "  -p   Primer fasta file"
	echo ""
	echo "Optional:"
	echo "  -e   Number of allowed mismatches near the 3' end, default: 1"
	echo "  -o   Output directory, default: ./virtualPCR_output_<yyyyMMddhhmmss_xxxxxx>"
	echo "  -k   Keep original report.out"
	echo "  -l   Amplicon length range: <minlen> <maxlen>, default: 50 1000"
	echo "  -m   Java memory in GB: <Xms> <Xmx>, default: 4 16"
	echo "  -h   Show this help"
	echo ""

	echo "Output files:"
	echo "  One summary file for all genomes in the output directory:"
	echo "  virtualPCR_summary_<e>bpMismatch_<min>-<max>.vPCR.out"
	echo "  If no amplicon output is present for a genome in report.out, summary writes 'no output'."
	echo "  Tested with virtualPCR.jar from commit 736e57c9282775e24ba4f145b04a670f31360fb8."
	echo ""
}

die() {
	echo "ERROR: $1" >&2
	exit 1
}

cleanup() {
	[[ -n "${run_config}" && -f "${run_config}" ]] && rm -f "${run_config}"
	[[ -n "${temp_input_dir}" && -d "${temp_input_dir}" ]] && rm -rf "${temp_input_dir}"
	[[ -n "${temp_java_output_path}" && -d "${temp_java_output_path}" ]] && rm -rf "${temp_java_output_path}"
	[[ -n "${temp_run_dir}" && -d "${temp_run_dir}" ]] && rmdir "${temp_run_dir}" 2>/dev/null || true
}

is_non_negative_int() {
	[[ "$1" =~ ^[0-9]+$ ]]
}

is_positive_int() {
	[[ "$1" =~ ^[1-9][0-9]*$ ]]
}

parse_args() {
	while [[ $# -gt 0 ]]; do
		case "$1" in
			-i)
				shift
				[[ $# -gt 0 ]] || die "-i requires at least one genome file"
				while [[ $# -gt 0 && "$1" != -* ]]; do
					input_files+=("$1")
					shift
				done
				continue
				;;
			-p)
				[[ $# -ge 2 ]] || die "-p requires a primer file"
				primer_file="$2"
				shift 2
				;;
			-e)
				[[ $# -ge 2 ]] || die "-e requires an integer value"
				number3error="$2"
				shift 2
				;;
			-o)
				[[ $# -ge 2 ]] || die "-o requires an output directory"
				output_root="$2"
				shift 2
				;;
			-k)
				keep_report=true
				shift
				;;
			-l)
				[[ $# -ge 3 ]] || die "-l requires two integers: minlen maxlen"
				minlen="$2"
				maxlen="$3"
				shift 3
				;;
			-m)
				[[ $# -ge 3 ]] || die "-m requires two integers: Xms Xmx"
				xms_gb="$2"
				xmx_gb="$3"
				shift 3
				;;
			-h|--help)
				show_help
				exit 0
				;;
			*)
				die "Unknown option: $1"
				;;
		esac
	done
}

main() {
	trap cleanup EXIT

	parse_args "$@"

	[[ ${#input_files[@]} -gt 0 ]] || {
		show_help
		die "At least one genome file is required via -i"
	}
	[[ -n "${primer_file}" ]] || {
		show_help
		die "Primer file is required via -p"
	}
	[[ -f "${JAR_PATH}" ]] || die "Missing jar file: ${JAR_PATH}"
	[[ -f "${primer_file}" ]] || die "Primer file does not exist: ${primer_file}"

	for genome in "${input_files[@]}"; do
		[[ -f "${genome}" ]] || die "Genome file does not exist: ${genome}"
	done

	is_non_negative_int "${number3error}" || die "-e must be a non-negative integer"

	is_non_negative_int "${minlen}" || die "-l minlen must be a non-negative integer"
	is_non_negative_int "${maxlen}" || die "-l maxlen must be a non-negative integer"
	(( minlen <= maxlen )) || die "-l requires minlen <= maxlen"

	is_positive_int "${xms_gb}" || die "-m Xms must be a positive integer"
	is_positive_int "${xmx_gb}" || die "-m Xmx must be a positive integer"
	(( xms_gb <= xmx_gb )) || die "-m requires Xms <= Xmx"

	current_time=$(date +"%Y%m%d%H%M%S")
	random_string=$(tr -dc 'a-z0-9' < /dev/urandom | head -c 6)

	if [[ -z "${output_root}" ]]; then
		output_root="./virtualPCR_output_${current_time}_${random_string}"
	fi

	mkdir -p "${output_root}" || die "Failed to create output directory: ${output_root}"
	output_root_abs="$(cd "${output_root}" && pwd)"

	temp_run_dir="${WORK_DIR}/virtualpcr_run_${current_time}_${random_string}"
	temp_input_dir="${temp_run_dir}/input"
	temp_java_output_dir="virtualpcr_java_output_${current_time}_${random_string}"
	temp_java_output_path="${WORK_DIR}/${temp_java_output_dir}"
	run_config="${temp_run_dir}/config.run.txt"
	run_report="${temp_java_output_path}/report.out"

	mkdir -p "${temp_input_dir}" || die "Failed to create temporary input directory"
	mkdir -p "${temp_java_output_path}" || die "Failed to create temporary Java output directory"
	cat > "${run_config}" << 'EOF' || die "Failed to create active run config"
# --- Paths ---
targets_path=
primers_path=
output_path=

# --- Search mode ---
type=primer
molecular=linear
linkedsearch=false
FRpairs=false
CTconversion=false

# --- Filters ---
minlen=50
maxlen=1000
number3errors=1

# --- Output control ---
primerstatistic=false
SequenceExtract=false
ShowPrimerAlignment=false
ShowPCRProducts=true
ShowOnlyAmplicons=false
ShowPrimerAlignmentPCRproduct=false
EOF

	declare -A amplicon_size

	# When target_path is set to a directory, virtualPCR.jar will process all files, including any unrelated files. 
	# To allow the screening of more than 1 genome file, all input file are sym linked to a temporary input directory, which is then set as the target_path in the run config.
	for genome in "${input_files[@]}"; do
		genome_base="$(basename "${genome}")"
		ln -s "$(realpath "${genome}")" "${temp_input_dir}/${genome_base}" || die "Failed to link ${genome_base}"
		amplicon_size["${genome_base}"]="no output"
	done

	sed -i -E "s|^targets_path=.*$|targets_path=${temp_input_dir}|" "${run_config}"
	sed -i -E "s|^primers_path=.*$|primers_path=$(realpath "${primer_file}")|" "${run_config}"
	sed -i -E "s|^output_path=.*$|output_path=${temp_java_output_dir}|" "${run_config}" # virtualPCR.jar output are forced lowercase. Saved to temporary directory to circumvent for case-sensitive file systems.
	sed -i -E "s|^number3errors=.*$|number3errors=${number3error}|" "${run_config}"
	sed -i -E "s|^minlen=.*$|minlen=${minlen}|" "${run_config}"
	sed -i -E "s|^maxlen=.*$|maxlen=${maxlen}|" "${run_config}"

	java_output=$(cd "${WORK_DIR}" && java -Xms"${xms_gb}"g -Xmx"${xmx_gb}"g -jar "${JAR_PATH}" "${run_config}" 2>&1)
	java_status=$?
	if [[ ${java_status} -ne 0 ]]; then
		echo "${java_output}" >&2
		die "Java run failed"
	fi

	[[ -f "${run_report}" ]] || die "Missing report.out"

	current_genome=""
	while IFS= read -r line; do
		if [[ "${line}" =~ ^In[[:space:]]silico[[:space:]]PCR,[[:space:]]primers[[:space:]]search[[:space:]]for:[[:space:]](.+)$ ]]; then
			current_path="${BASH_REMATCH[1]}"
			current_genome="$(basename "${current_path}")"
			continue
		fi

		if [[ -n "${current_genome}" && -n "${amplicon_size[${current_genome}]+x}" ]]; then
			size_bp=""
			if [[ "${line}" =~ ^Amplicon[[:space:]]size:[[:space:]]*([0-9]+)bp$ ]]; then
				size_bp="${BASH_REMATCH[1]}bp"
			fi

			if [[ -z "${size_bp}" ]]; then
				continue
			fi

			if [[ "${amplicon_size[${current_genome}]}" == "no output" ]]; then
				amplicon_size["${current_genome}"]="${size_bp}"
			else
				amplicon_size["${current_genome}"]+=",${size_bp}"
			fi
		fi
	done < "${run_report}"

	summary_out="${output_root_abs}/virtualPCR_summary_${number3error}bpMismatch_${minlen}-${maxlen}.vPCR.out"
	{
		echo -e "Genome\tAmpliconSizes"
	for genome in "${input_files[@]}"; do
		genome_base="$(basename "${genome}")"
		echo -e "${genome_base}\t${amplicon_size[${genome_base}]}"
	done
	} > "${summary_out}"
	echo "Saved: ${summary_out}"

	if [[ "${keep_report}" == true ]]; then
		mv "${run_report}" "${output_root_abs}/report.out" || die "Failed to move report.out"
		echo "Saved: ${output_root_abs}/report.out"
	fi

	shopt -s dotglob nullglob
	for java_item in "${temp_java_output_path}"/*; do
		item_name="$(basename "${java_item}")"
		if [[ "${item_name}" == "report.out" ]]; then
			continue
		fi
		mv "${java_item}" "${output_root_abs}/" || die "Failed to move Java output item: ${item_name}"
	done
	shopt -u dotglob nullglob

	if [[ "${keep_report}" != true ]]; then
		rm -f "${run_report}"
	fi
}

input_files=()
primer_file=""
number3error="1"
output_root=""
keep_report=false
minlen="50"
maxlen="1000"
xms_gb="4"
xmx_gb="16"

main "$@"
