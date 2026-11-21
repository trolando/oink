# Copyright 2026 Tom van Dijk
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Smoke tests of bench_solvers: short headless runs on a few small games,
# checking the exit status, the reports and the output. Run as
#
#   cmake -DBENCH=<bench_solvers> -DDATA=<tests dir> -DWORK=<dir> -DMODE=<mode> -P smoke.cmake
#
# where MODE selects one of the runs below.

if(NOT BENCH OR NOT DATA OR NOT WORK OR NOT MODE)
    message(FATAL_ERROR "usage: cmake -DBENCH=... -DDATA=... -DWORK=... -DMODE=... -P smoke.cmake")
endif()

file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
string(ASCII 27 ESC)
set(GAMES "${DATA}/vb001" "${DATA}/vb002" "${DATA}/vb003")

# Run bench_solvers with the remaining arguments in WORK; check the exit status.
function(bench expected)
    execute_process(COMMAND "${BENCH}" ${ARGN}
        WORKING_DIRECTORY "${WORK}"
        RESULT_VARIABLE status
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err)
    if(NOT status STREQUAL "${expected}")
        message(FATAL_ERROR "bench_solvers ${ARGN}\nexit status ${status}, expected ${expected}\nstdout:\n${out}\nstderr:\n${err}")
    endif()
    set(out "${out}" PARENT_SCOPE)
    set(err "${err}" PARENT_SCOPE)
endfunction()

function(expect text haystack what)
    string(FIND "${haystack}" "${text}" at)
    if(at EQUAL -1)
        message(FATAL_ERROR "${what}: missing \"${text}\" in:\n${haystack}")
    endif()
endfunction()

function(reject text haystack what)
    string(FIND "${haystack}" "${text}" at)
    if(NOT at EQUAL -1)
        message(FATAL_ERROR "${what}: unexpected \"${text}\" in:\n${haystack}")
    endif()
endfunction()

# Count the matches of <regex> in <text> into <var>.
function(count var regex text)
    string(REGEX MATCHALL "${regex}" matches "${text}")
    list(LENGTH matches n)
    set(${var} ${n} PARENT_SCOPE)
endfunction()

# The sequence of trials of a JSON Lines report: game, solver and phase.
function(schedule var file)
    file(READ "${file}" text)
    string(REGEX MATCHALL "\"game_index\":[0-9]+,\"solver\":\"[^\"]*\",\"phase\":\"[a-z]+\"" trials "${text}")
    set(${var} "${trials}" PARENT_SCOPE)
endfunction()

if(MODE STREQUAL "basic")
    # Two solvers, three games, one warm-up and two measured rounds.
    bench(0 --tl --zlk --repeat=2 --warmup=1 --min-solve-time=0 --bench-seed=7 --jsonl=run.jsonl --csv=run.csv ${GAMES})
    reject("${ESC}" "${out}${err}" "no escape sequences without a terminal")
    expect("Benchmark summary" "${out}" "summary")
    expect("VBS W/T" "${out}" "summary")
    expect("Schedule seed: 7" "${out}" "summary")
    reject("geom. speedup" "${out}" "no comparison without a reference")
    expect("game\ttl\tzlk" "${err}" "per-game header")
    file(STRINGS "${WORK}/run.jsonl" lines)
    list(GET lines 0 first)
    list(GET lines -1 last)
    expect("\"event\":\"metadata\"" "${first}" "metadata first")
    expect("\"event\":\"end\"" "${last}" "end marker last")
    file(READ "${WORK}/run.jsonl" jsonl)
    count(samples "\"event\":\"sample\"" "${jsonl}")
    count(ok "\"status\":\"ok\"" "${jsonl}")
    if(NOT samples EQUAL 18)
        message(FATAL_ERROR "expected 18 samples (3 games, 2 solvers, 3 rounds), found ${samples}")
    endif()
    count(results "\"event\":\"result\"" "${jsonl}")
    if(NOT results EQUAL 6)
        message(FATAL_ERROR "expected 6 results, found ${results}")
    endif()
    file(STRINGS "${WORK}/run.csv" rows)
    list(LENGTH rows n)
    if(NOT n EQUAL 7)
        message(FATAL_ERROR "expected a CSV header and 6 rows, found ${n} lines")
    endif()
    list(GET rows 0 header)
    expect("game,graph_fnv1a64,vertices,edges" "${header}" "CSV header")

    # The reports are not replaced without --overwrite.
    bench(3 --tl --quick --jsonl=run.jsonl ${GAMES})
    expect("already exists" "${err}" "existing report")

    # The same schedule seed gives the same schedule.
    bench(0 --tl --zlk --repeat=2 --warmup=1 --min-solve-time=0 --bench-seed=7 --jsonl=again.jsonl ${GAMES})
    schedule(a "${WORK}/run.jsonl")
    schedule(b "${WORK}/again.jsonl")
    if(NOT a STREQUAL b)
        message(FATAL_ERROR "the same --bench-seed gave another schedule")
    endif()

elseif(MODE STREQUAL "sampling")
    # Four solves under the default deadline cannot reach this target.
    bench(0 --tl --repeat=2 --warmup=0 --min-solve-time=200 --max-repeat=4 --bench-seed=3 --jsonl=run.jsonl --csv=run.csv "${DATA}/vb001")
    file(READ "${WORK}/run.jsonl" jsonl)
    count(measured "\"phase\":\"measured\"" "${jsonl}")
    if(NOT measured EQUAL 4)
        message(FATAL_ERROR "expected 4 measured trials, found ${measured}")
    endif()
    expect("\"status\":\"ok\",\"failure_phase\":null,\"samples\":4" "${jsonl}" "four verified solves")
    expect("\"sampling_capped\":true,\"mean_seconds\":" "${jsonl}" "cap and mean")
    file(READ "${WORK}/run.csv" csv)
    expect("mean_seconds,accumulated_solve_seconds,sampling_capped" "${csv}" "CSV sampling fields")
    expect(",true," "${csv}" "CSV cap")
    expect("n=4" "${err}" "displayed repetition count")

elseif(MODE STREQUAL "reference")
    bench(0 --tl --tlq --zlk --quick --reference=tl --bench-seed=1 --jsonl=run.jsonl ${GAMES})
    expect("Reference: tl." "${out}" "summary")
    expect("geom. speedup" "${out}" "comparison columns")
    expect("W/T/L" "${out}" "comparison columns")
    file(READ "${WORK}/run.jsonl" jsonl)
    expect("\"reference\":\"tl\"" "${jsonl}" "reference in the report")
    # An unselected or ambiguous reference is an error.
    bench(3 --tl --quick --reference=zlk ${GAMES})
    expect("not selected" "${err}" "unselected reference")

elseif(MODE STREQUAL "random")
    bench(0 --tl --zlk --fpi --psi --tlpar=parmin=1,chunk=4 --workers=2 --quick --count=3 --size=20 --bench-seed=11 --jsonl=run.jsonl)
    count(games "seed=" "${err}")
    if(NOT games EQUAL 3)
        message(FATAL_ERROR "expected 3 random games:\n${err}")
    endif()
    expect("Game-series seed:" "${out}" "series seed")
    file(READ "${WORK}/run.jsonl" jsonl)
    foreach(id zlk fpi psi tlpar)
        count(parallel "\"parallel\":true,\"effective_solver\":\"${id}\"" "${jsonl}")
        if(NOT parallel EQUAL 3)
            message(FATAL_ERROR "expected ${id} to run in parallel on all 3 small games")
        endif()
    endforeach()
    expect("\"effective_solver\":\"tlpar\",\"effective_configuration\":\"parmin=1,chunk=4\"" "${jsonl}" "parallel TLPar configuration")

elseif(MODE STREQUAL "errors")
    bench(3 --tlpar ${GAMES})
    expect("solver tlpar is parallel-only" "${err}" "explicit sequential parallel-only solver")
    bench(3 --all --tlpar --force-parallel ${GAMES})
    expect("solver tlpar is parallel-only" "${err}" "explicit parallel-only choice with all and no workers")
    bench(0 --quick --jsonl=default.jsonl "${DATA}/vb001")
    expect("Sequential run omits parallel-only solver tlpar" "${err}" "implicit sequential exclusion")
    file(READ "${WORK}/default.jsonl" defaults)
    reject("\"solver\":\"tlpar\"" "${defaults}" "no parallel-only solver in sequential defaults")
    bench(3 --tlpar=bogus --workers=2 ${GAMES})
    expect("invalid configuration" "${err}" "invalid configuration")
    bench(3 --no-such-option ${GAMES})
    expect("unknown option" "${err}" "unknown option")
    bench(3 --tl "${WORK}/missing-file")
    expect("not found" "${err}" "missing input")
    # an unreadable game is a failure of every solver, not the end of the run
    file(WRITE "${WORK}/broken.pg" "parity 3;\n0 0 0 9;\n")
    bench(1 --tl --zlk --quick --sort "${WORK}/broken.pg" "${DATA}/vb001")
    expect("input_error" "${err}" "unreadable game")

elseif(MODE STREQUAL "preprocess")
    # once-per-game preprocessing; a single-parity game needs no solver at all
    file(WRITE "${WORK}/even.pg" "parity 2;\n0 2 0 1;\n1 4 1 0,2;\n2 0 1 2;\n")
    # vb003, vb004 and vb006 keep a residual; vb001 and vb002 are solved by preprocessing
    bench(0 --preprocess --tl --fpjm --quick --sort --jsonl=run.jsonl --csv=run.csv "${DATA}/vb003" "${DATA}/vb004" "${DATA}/vb006" "${WORK}/even.pg")
    expect("Preprocessing: ON" "${out}" "summary")
    expect("Fully solved by preprocessing: 1/4" "${out}" "summary")
    expect("preprocessed" "${err}" "per-game line")
    file(READ "${WORK}/run.jsonl" jsonl)
    count(preparations "\"event\":\"preparation\",\"game_index\":[0-9]+,\"enabled\":true,\"status\":\"ok\"" "${jsonl}")
    if(NOT preparations EQUAL 4)
        message(FATAL_ERROR "expected 4 preparations, found ${preparations}")
    endif()
    count(preprocessed "\"status\":\"preprocessed\"" "${jsonl}")
    if(NOT preprocessed EQUAL 2)
        message(FATAL_ERROR "expected the results of 2 solvers to be preprocessed, found ${preprocessed}")
    endif()
    reject("\"status\":\"invalid_solution\"" "${jsonl}" "lifted solutions verify")
    bench(3 --no-preprocess --single --tl ${GAMES})
    expect("preprocessing is off" "${err}" "contradicting options")

elseif(MODE STREQUAL "live")
    # without a terminal the live table falls back to plain output
    bench(0 --live --tl --quick ${GAMES})
    reject("${ESC}" "${out}${err}" "plain output without a terminal")
    expect("game\ttl" "${err}" "plain per-game lines")
    bench(3 --review --tl ${GAMES})
    expect("--review needs" "${err}" "review without the live table")

elseif(MODE STREQUAL "experiment")
    # a named experiment saves its plan, games and reports
    bench(0 --experiment=exp --tl --zlk --quick --sort --bench-seed=5 ${GAMES})
    file(READ "${WORK}/exp/status.json" status)
    expect("\"complete\":true" "${status}" "status")
    file(GLOB blocks "${WORK}/exp/games/*.jsonl")
    list(LENGTH blocks n)
    if(NOT n EQUAL 3)
        message(FATAL_ERROR "expected 3 finished games, found ${n}")
    endif()
    file(STRINGS "${WORK}/exp/results.csv" rows)
    list(LENGTH rows n)
    if(NOT n EQUAL 7)
        message(FATAL_ERROR "expected a CSV header and 6 rows, found ${n} lines")
    endif()
    file(READ "${WORK}/exp/results.jsonl" results)
    count(games "\"event\":\"game\"" "${results}")
    count(summaries "\"event\":\"summary\"" "${results}")
    if(NOT games EQUAL 3 OR NOT summaries EQUAL 2)
        message(FATAL_ERROR "results.jsonl: ${games} games and ${summaries} summaries")
    endif()

    # resuming a finished experiment measures nothing
    bench(0 --resume=exp)
    file(READ "${WORK}/exp/sessions/2.jsonl" session)
    reject("\"event\":\"sample\"" "${session}" "no new trials")
    expect("(3/3 games)" "${out}" "experiment results")

    # an existing experiment is not started again; --resume takes no other options
    bench(3 --experiment=exp --tl ${GAMES})
    expect("--resume" "${err}" "existing experiment without a terminal")
    bench(3 --resume=exp --tl)
    expect("only --plain" "${err}" "options with --resume")
    bench(3 --experiment=other --jsonl=x.jsonl --tl ${GAMES})
    expect("writes its own reports" "${err}" "reports of an experiment")

    # a game that is not in done.txt is measured again; a changed input is refused
    file(MAKE_DIRECTORY "${WORK}/in")
    file(COPY "${DATA}/vb001" "${DATA}/vb002" DESTINATION "${WORK}/in")
    bench(0 --experiment=exp2 --tl --quick --sort "${WORK}/in")
    file(STRINGS "${WORK}/exp2/done.txt" finished)
    list(GET finished 0 first)
    file(WRITE "${WORK}/exp2/done.txt" "${first}\ngame 1 \"truncated")
    bench(0 --resume=exp2)
    file(READ "${WORK}/exp2/sessions/2.jsonl" session)
    count(trials "\"event\":\"sample\"" "${session}")
    if(NOT trials EQUAL 1)
        message(FATAL_ERROR "expected the one unfinished game to run again, found ${trials} trials")
    endif()
    bench(0 --resume=exp2)
    file(READ "${WORK}/exp2/sessions/3.jsonl" session)
    reject("\"event\":\"sample\"" "${session}" "the new commit survived a truncated record")
    file(WRITE "${WORK}/exp2/done.txt" "${first}\n")
    file(APPEND "${WORK}/in/vb002" "\n")
    bench(3 --resume=exp2)
    expect("an input changed" "${err}" "changed input")

elseif(MODE STREQUAL "timeout")
    # POSIX: an external solver that hangs is stopped by the timeout, and a
    # crash of the worker is caught; the other solvers continue.
    bench(1 --tl -e "slow:sleep 5" -e "boom:kill -9 $PPID" --repeat=1 --warmup=1 --min-solve-time=0 --timeout=0.3 --bench-seed=2 --jsonl=run.jsonl --csv=run.csv --artifacts=failed "${DATA}/vb001" "${DATA}/vb002")
    file(READ "${WORK}/run.jsonl" jsonl)
    expect("\"status\":\"timeout\",\"failure_phase\":\"solve\"" "${jsonl}" "solve timeout")
    expect("\"status\":\"crash\",\"failure_phase\":\"solve\"" "${jsonl}" "solve crash")
    file(STRINGS "${WORK}/run.csv" timeouts REGEX ",timeout,.*[,]solve,,,[0-9]*,1,0,")
    list(LENGTH timeouts timeout_count)
    if(NOT timeout_count EQUAL 2)
        message(FATAL_ERROR "expected two CSV rows with a solve timeout")
    endif()
    expect("\"status\":\"skipped\"" "${jsonl}" "later trials skipped")
    count(ok "\"solver\":\"tl\",\"phase\":\"[a-z]+\",[^\n]*\"status\":\"ok\"" "${jsonl}")
    if(NOT ok EQUAL 4)
        message(FATAL_ERROR "tl should succeed in all 4 trials, found ${ok}")
    endif()
    file(GLOB_RECURSE saved "${WORK}/failed/*.pg")
    if(NOT saved)
        message(FATAL_ERROR "no failure artifacts saved")
    endif()

elseif(MODE STREQUAL "skip")
    file(WRITE "${WORK}/flaky.sh" [=[
case $(cat flaky-count 2>/dev/null) in
    '') echo 1 > flaky-count; sleep 5 ;;
    1) echo 2 > flaky-count; sleep 5 ;;
    *) "$1" "$2" "$3" ;;
esac
]=])
    foreach(mode reuse fresh)
        set(isolation)
        if(mode STREQUAL "fresh")
            set(isolation --fresh-process)
        endif()
        # The warm-up times out too. Each pair still gets two measured attempts.
        bench(1 ${isolation} --tl -e "slow:sleep 5" --repeat=4 --warmup=1
            --min-solve-time=0 --timeout=0.5 --bench-seed=7 --artifacts=
            --jsonl=${mode}.jsonl --csv=${mode}.csv "${DATA}/vb001" "${DATA}/vb002")
        file(READ "${WORK}/${mode}.jsonl" jsonl)
        count(real "\"solver\":\"slow\",\"phase\":\"measured\"[^\n]*\"status\":\"timeout\"[^\n]*\"executed\":true" "${jsonl}")
        count(skipped "\"solver\":\"slow\",\"phase\":\"measured\"[^\n]*\"status\":\"timeout\"[^\n]*\"executed\":false,\"skip_reason\":\"timeout_limit\"[^\n]*not repeated after two timeouts" "${jsonl}")
        count(fast "\"solver\":\"tl\",\"phase\":\"measured\"[^\n]*\"status\":\"ok\"" "${jsonl}")
        if(NOT real EQUAL 4 OR NOT skipped EQUAL 4 OR NOT fast EQUAL 8)
            message(FATAL_ERROR "${mode}: ${real} real timeouts, ${skipped} skipped timeouts, ${fast} fast trials")
        endif()
        file(STRINGS "${WORK}/${mode}.csv" scores REGEX ",timeout,0,.*[,]2,2,1$")
        list(LENGTH scores n)
        if(NOT n EQUAL 2)
            message(FATAL_ERROR "${mode}: expected two timed-out pairs with PAR2 of 1 second")
        endif()

        # The first two invocations (one warm-up and one measured) time out.
        file(REMOVE "${WORK}/flaky-count")
        bench(1 ${isolation} --tl -e "flaky:sh flaky.sh \"${OINK}\" %I %O" --repeat=4 --warmup=1
            --min-solve-time=0 --timeout=0.5 --bench-seed=7 --artifacts=
            --jsonl=${mode}-flaky.jsonl "${DATA}/vb001")
        file(READ "${WORK}/${mode}-flaky.jsonl" flaky)
        expect("\"samples\":3,\"target_samples\":4,\"solve_timeouts\":1,\"skipped_timeouts\":0" "${flaky}" "retries after one measured timeout")
    endforeach()

    bench(1 -e "slow:sleep 5" --repeat=3 --warmup=0 --min-solve-time=0 --timeout=0.2
        --skip-after-timeouts=0 --artifacts= --jsonl=all.jsonl "${DATA}/vb001")
    file(READ "${WORK}/all.jsonl" all)
    expect("\"solve_timeouts\":3,\"skipped_timeouts\":0" "${all}" "skip disabled")
    bench(3 --tl --skip-after-timeouts=-1 "${DATA}/vb001")

elseif(MODE STREQUAL "fresh")
    set(args --tl --zlk --repeat=2 --warmup=1 --min-solve-time=0 --bench-seed=7)
    bench(0 ${args} --jsonl=reuse.jsonl --csv=reuse.csv ${GAMES})
    bench(0 --fresh-process ${args} --jsonl=fresh.jsonl --csv=fresh.csv ${GAMES})
    schedule(reuse "${WORK}/reuse.jsonl")
    schedule(fresh "${WORK}/fresh.jsonl")
    if(NOT reuse STREQUAL fresh)
        message(FATAL_ERROR "fresh processes changed the schedule")
    endif()
    file(READ "${WORK}/reuse.jsonl" reuse)
    expect("\"isolation\":\"worker process per game, restarted after a failure\"" "${reuse}" "reuse is the default")
    reject("\"worker_generation\":2" "${reuse}" "one worker per game by default")
    reject("\"trial_peak_rss_bytes\":" "${reuse}" "no per-trial RSS in reuse mode")
    file(READ "${WORK}/fresh.jsonl" fresh)
    expect("\"isolation\":\"new executable per trial\"" "${fresh}" "fresh mode metadata")
    file(STRINGS "${WORK}/fresh.jsonl" samples REGEX "\"event\":\"sample\"")
    set(expected_game 0)
    set(generation 1)
    foreach(sample IN LISTS samples)
        expect("\"game_index\":${expected_game}," "${sample}" "fresh sample game")
        expect("\"worker_generation\":${generation}," "${sample}" "new process for each trial")
        expect("\"status\":\"ok\"" "${sample}" "verified trial")
        if(NOT sample MATCHES "\"load_seconds\":[0-9][^,]*,\"parse_sort_seconds\":[0-9]")
            message(FATAL_ERROR "missing load or parse time: ${sample}")
        endif()
        if(CMAKE_HOST_APPLE OR CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
            if(NOT sample MATCHES "\"trial_peak_rss_bytes\":([1-9][0-9]*),\"solve_start_peak_rss_bytes\":([1-9][0-9]*),\"solve_end_peak_rss_bytes\":([1-9][0-9]*)")
                message(FATAL_ERROR "missing fresh RSS readings: ${sample}")
            endif()
            if(CMAKE_MATCH_1 LESS CMAKE_MATCH_3 OR CMAKE_MATCH_3 LESS CMAKE_MATCH_2)
                message(FATAL_ERROR "process high-water mark went down: ${sample}")
            endif()
        endif()
        math(EXPR generation "${generation}+1")
        if(generation EQUAL 7)
            set(generation 1)
            math(EXPR expected_game "${expected_game}+1")
        endif()
    endforeach()
    if(NOT expected_game EQUAL 3 OR NOT generation EQUAL 1)
        message(FATAL_ERROR "expected six trials for each of three games")
    endif()
    if(CMAKE_HOST_APPLE OR CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        file(STRINGS "${WORK}/reuse.csv" reuse_rss REGEX ",,,[1-9][0-9]*,0,0,")
        file(STRINGS "${WORK}/fresh.csv" fresh_rss REGEX ",[1-9][0-9]*,[1-9][0-9]*,,0,0,")
        list(LENGTH reuse_rss reuse_rows)
        list(LENGTH fresh_rss fresh_rows)
        if(NOT reuse_rows EQUAL 6 OR NOT fresh_rows EQUAL 6)
            message(FATAL_ERROR "CSV RSS columns did not match the isolation mode")
        endif()
    endif()

    bench(0 --fresh-process --workers=2 --tlpar=parmin=1,chunk=4 --quick --count=1 --size=20 --seed=23 --jsonl=random.jsonl)
    file(READ "${WORK}/random.jsonl" random)
    expect("\"load_seconds\":null,\"parse_sort_seconds\":" "${random}" "generated input has no file load")
    expect("\"parallel\":true,\"effective_solver\":\"tlpar\",\"effective_configuration\":\"parmin=1,chunk=4\",\"actual_workers\":2" "${random}" "fresh Lace workers and configuration")
    bench(0 --fresh-process --experiment=exp --preprocess --tl --quick "${DATA}/vb003")
    file(READ "${WORK}/exp/results.csv" saved_csv)
    bench(0 --resume=exp)
    file(READ "${WORK}/exp/results.csv" resumed_csv)
    if(NOT saved_csv STREQUAL resumed_csv)
        message(FATAL_ERROR "fresh-process statistics changed on resume")
    endif()
    file(READ "${WORK}/exp/results.jsonl" resumed)
    expect("\"isolation\":\"new executable per trial\"" "${resumed}" "fresh mode survives resume")
    expect("\"once_per_game\":false" "${resumed}" "fresh preprocessing schedule")
    expect("\"included_in_trial_setup\":true" "${resumed}" "fresh preprocessing is measured in setup")

    # Simulate a Lace allocation failure followed by a fatal worker exit.
    bench(1 --fresh-process --tl -e "slow:sleep 5"
        -e "memory:printf 'Lace error: unable to allocate memory for the workers\\n' >&2 && kill -9 $PPID"
        --repeat=2 --warmup=0 --min-solve-time=0 --timeout=0.3 --bench-seed=2 --jsonl=failures.jsonl "${DATA}/vb001")
    file(READ "${WORK}/failures.jsonl" failed)
    expect("\"status\":\"timeout\",\"failure_phase\":\"solve\"" "${failed}" "fresh solve deadline")
    expect("\"status\":\"memory_failure\",\"failure_phase\":\"solve\"" "${failed}" "fresh memory failure")
    expect("\"status\":\"skipped\"" "${failed}" "fresh failed cell skips later trials")
    count(ok "\"solver\":\"tl\",\"phase\":\"measured\",[^\n]*\"status\":\"ok\"" "${failed}")
    if(NOT ok EQUAL 2)
        message(FATAL_ERROR "other solvers should finish after a fresh worker fails")
    endif()
    bench(1 --fresh-process --tl --quick --safety-timeout=0.000001 --jsonl=load.jsonl "${DATA}/vb001")
    file(READ "${WORK}/load.jsonl" load)
    expect("\"status\":\"safety_timeout\",\"failure_phase\":\"load\"" "${load}" "fresh load deadline")
    file(WRITE "${WORK}/broken.pg" "parity 3;\n0 0 0 9;\n")
    bench(1 --fresh-process --tl --quick --jsonl=broken.jsonl "${WORK}/broken.pg")
    file(READ "${WORK}/broken.jsonl" broken)
    expect("\"status\":\"input_error\",\"failure_phase\":\"load\"" "${broken}" "fresh parser failure")
    bench(3 --fresh-process --in-process --tl --quick "${DATA}/vb001")
    expect("choose --fresh-process or --in-process" "${err}" "incompatible process modes")

else()
    message(FATAL_ERROR "unknown MODE ${MODE}")
endif()
