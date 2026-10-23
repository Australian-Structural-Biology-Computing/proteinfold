#!/usr/bin/env bash

original_args=("$@")
command_name="${1:-}"
image=''
cache_before=''

if [[ "$command_name" == run || "$command_name" == create ]]; then
    shift
    while (( $# )); do
        case "$1" in
            --env|--env-file|--volume|--mount|--workdir|--user|--name|--network|--entrypoint|--hostname|--platform|--pull|--memory|--cpus|--gpus|--device|--label|--publish|--tmpfs|--add-host|--group-add|--ulimit|--cap-add|--cap-drop|--security-opt|--shm-size|--ipc|--pid|--cpu-shares|--cpuset|--cpuset-cpus|--memory-swap|--oom-score-adj|--cgroup-parent|--sysctl|--dns|--dns-search|--dns-option|--log-driver|--log-opt|--restart|--stop-signal|--health-cmd|--volumes-from|--link|-e|-v|-w|-u|-p|-m|-l|-h)
                shift 2
                ;;
            --*=*|-e?*|-v?*|-w?*|-u?*|-p?*|-m?*|-l?*|-h?*)
                shift
                ;;
            -* )
                shift
                ;;
            *)
                if [[ "$1" == */* || "$1" == *:* ]]; then
                    image="$1"
                fi
                break
                ;;
        esac
    done
    if [[ -n "$image" ]]; then
        if "$DOCKER_CI_BINARY" image inspect "$image" >/dev/null 2>&1; then
            cache_before=present
        else
            cache_before=missing
        fi
    fi
elif [[ "$command_name" == pull ]]; then
    image="${*: -1}"
fi

started=$(date +%s%N)
"$DOCKER_CI_BINARY" "${original_args[@]}"
status=$?
ended=$(date +%s%N)

case "$command_name" in
    pull|create|start|run)
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$command_name" "$(( (ended - started) / 1000000 ))" "$status" "$image" "$cache_before" >> "$DOCKER_CI_TIMING_LOG"
        ;;
esac

exit "$status"
