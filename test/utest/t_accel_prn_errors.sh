#!/usr/bin/env bash
# t_accel_prn_errors.sh : exercise accel_prn FATAL paths via rnx2rtkp
# Each case feeds a deliberately bad CSV and asserts that rnx2rtkp:
#   (a) exits non-zero
#   (b) prints the expected "accel_prn: FATAL: ...:N <reason>" message
#
# Tests covered:
#   A1  empty file
#   A2  header-only (no data rows)
#   A5  duplicate timestamp
#   A6  reverse order
#   A8  unparseable time
# (A3/A4/A7 dropped: trailing columns are intentionally tolerated, and the
#  CSV no longer carries prn values so there's no "negative" path.)
set -u

BIN="${RNX2RTKP:-../../app/consapp/rnx2rtkp/gcc/rnx2rtkp}"
[[ -x "$BIN" ]] || { echo "rnx2rtkp not found at $BIN; build it first" >&2; exit 2; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail_count=0
pass_count=0

assert_fatal() {
    local name="$1" csv="$2" expected_substr="$3"
    local out
    out="$("$BIN" -pn "$csv" /dev/null 2>&1)"
    local rc=$?
    if [[ $rc -eq 0 ]]; then
        echo "  FAIL [$name]: expected non-zero exit, got 0"
        fail_count=$((fail_count+1)); return
    fi
    if ! grep -qF -- "$expected_substr" <<<"$out"; then
        echo "  FAIL [$name]: missing expected substring '$expected_substr'"
        echo "       got: $out"
        fail_count=$((fail_count+1)); return
    fi
    echo "  ok   [$name]"
    pass_count=$((pass_count+1))
}

# A1: empty file
: > "$TMP/empty.csv"
assert_fatal "A1 empty"        "$TMP/empty.csv"        "empty file"

# A2: header only
echo "time" > "$TMP/header_only.csv"
assert_fatal "A2 header-only"  "$TMP/header_only.csv"  "no data rows"

# A5: duplicate timestamp (== violates strict-increasing)
cat > "$TMP/dup.csv" <<EOF
time
2024-01-15T12:00:00.000Z
2024-01-15T12:00:00.000Z
EOF
assert_fatal "A5 duplicate-time" "$TMP/dup.csv"        "not strictly increasing"

# A6: reverse order
cat > "$TMP/reverse.csv" <<EOF
time
2024-01-15T12:00:01.000Z
2024-01-15T12:00:00.000Z
EOF
assert_fatal "A6 reverse"      "$TMP/reverse.csv"      "not strictly increasing"

# A8: unparseable time
cat > "$TMP/badtime.csv" <<EOF
time
not-a-time
EOF
assert_fatal "A8 bad-time"     "$TMP/badtime.csv"      "unparseable time"

# A_extra: missing file
assert_fatal "A-extra missing-file" "$TMP/does_not_exist.csv" "cannot open"

echo
echo "passed: $pass_count, failed: $fail_count"
[[ $fail_count -eq 0 ]] || exit 1
exit 0
