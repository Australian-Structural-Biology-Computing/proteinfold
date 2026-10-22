#!/usr/bin/env bash

started=$(date +%s%N)
"$DOCKER_CI_BINARY" "$@"
status=$?
ended=$(date +%s%N)

case "${1:-}" in
    pull|create|start|run)
        image=''
        if [[ "$1" == pull ]]; then
            image="${*: -1}"
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$(( (ended - started) / 1000000 ))" "$status" "$image" >> "$DOCKER_CI_TIMING_LOG"
        ;;
esac

exit "$status"
