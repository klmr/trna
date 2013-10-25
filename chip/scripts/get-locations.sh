#/usr/bin/env bash

set -e
set -u

input=(results/de-genes/*.tsv)
reference='../common/data/tRNA.tsv'
prefix='de-genes-'

tsv-add-extension() {
    echo "$1"
}

tsv-write-header() {
    head -n 1 "$2" | tr '\n' '	' >> "$1"
    echo 'start	end' >> "$1"
}

tsv-write-body() {
    echo "$2" >> "$1"
}

bed-add-extension() {
    echo "${1%tsv}bed"
}

bed-write-header() {
    echo "track name=\"$(basename "${1%.bed}")\"" > "$1"
}

bed-write-body() {
    if [ "$2" == "" ]; then
        return
    fi
    while read no id bm bmA bmB fC l2fC p pa start end; do
        # Format of $id: "chrx.trna1234"
        id=$(tr -d '"' <<< $id)
        chr=$(cut -f 1 -d. <<< $id)
        trna=$(cut -f 2 -d. <<< $id)
        isoform=$(grep "$chr	 .$trna" "$reference" | cut -f 3)
        if [ $start -gt $end ]; then
            tmp=$start
            start=$end
            end=$tmp
            strand='-'
        else
            strand='+'
        fi
        # Simply scale p-value up to 1000; the higher the more significant
        score=$(perl -e "print int(1000 - $p * 20000)")
        echo "$chr	$start	$end	$isoform	$score	$strand"
    done <<< "$2" >> "$1"
}


case ${1---tsv} in
    --tsv)
        method=tsv
        ;;
    --bed)
        method=bed
        ;;
    *)
        echo >&2 "Invalid argument. Usage:"
        echo >&2
        echo >&2 "  $0 [--csv|--bed]"
        exit 1
esac

output="results/de-genes/$method"
mkdir -p "$output"

for file in "${input[@]}"; do
    base="$($method-add-extension "$(basename "$file")")"
    echo $base

    target="$output/$prefix$base"
    ids="$(tail -n+2 "$file" | cut -f 2 | sed 's/"//g')"
    rm -f "$target"
    $method-write-header "$target" "$file"
    contents=$(paste <(tail -n+2 "$file") <(for id in $ids; do
        id="$(sed 's/\./	 ./' <<< $id)	"
        grep "$id" "$reference" | cut -f 4,5
    done))
    $method-write-body "$target" "$contents"
done
