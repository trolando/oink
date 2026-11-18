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

elseif(MODE STREQUAL "calibrated")
    # A sampling target far above the solve time: three pilots, then more samples.
    bench(0 --tl --repeat=2 --warmup=0 --min-solve-time=0.0005 --max-repeat=8 --bench-seed=3 --jsonl=run.jsonl "${DATA}/vb001")
    file(READ "${WORK}/run.jsonl" jsonl)
    count(pilots "\"phase\":\"calibration\"" "${jsonl}")
    count(plans "\"event\":\"calibration\"" "${jsonl}")
    count(measured "\"phase\":\"measured\"" "${jsonl}")
    if(NOT pilots EQUAL 3 OR NOT plans EQUAL 1)
        message(FATAL_ERROR "expected 3 pilot trials and one calibration, found ${pilots} and ${plans}")
    endif()
    if(measured LESS 2 OR measured GREATER 8)
        message(FATAL_ERROR "expected 2 to 8 measured trials, found ${measured}")
    endif()

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
    bench(0 --tl --zlk --quick --count=3 --size=20 --bench-seed=11 --jsonl=run.jsonl)
    count(games "seed=" "${err}")
    if(NOT games EQUAL 3)
        message(FATAL_ERROR "expected 3 random games:\n${err}")
    endif()
    expect("Game-series seed:" "${out}" "series seed")

elseif(MODE STREQUAL "errors")
    bench(3 --tlpar=bogus ${GAMES})
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
    file(WRITE "${WORK}/exp2/done.txt" "${first}\n")
    bench(0 --resume=exp2)
    file(READ "${WORK}/exp2/sessions/2.jsonl" session)
    count(trials "\"event\":\"sample\"" "${session}")
    if(NOT trials EQUAL 1)
        message(FATAL_ERROR "expected the one unfinished game to run again, found ${trials} trials")
    endif()
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
    file(STRINGS "${WORK}/run.csv" timeouts REGEX ",timeout,.*[,]solve$")
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

else()
    message(FATAL_ERROR "unknown MODE ${MODE}")
endif()
