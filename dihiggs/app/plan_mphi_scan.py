#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import time
from pathlib import Path
from typing import Dict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Planificador de barridos en m_phi para PhysParamScan (2HDM)."
    )
    parser.add_argument(
        "--tan-beta",
        type=float,
        required=True,
        help="Valor de tan(beta).",
    )
    parser.add_argument(
        "--N-mphi",
        type=int,
        default=1,
        help="Parámetro N_mphi para PhysParamScan (por defecto 1).",
    )
    parser.add_argument(
        "--N-m12",
        type=int,
        required=True,
        help="Parámetro N_m12 para PhysParamScan.",
    )
    parser.add_argument(
        "--tol-exp",
        type=int,
        default=None,
        help="tol_exp (opcional).",
    )
    parser.add_argument(
        "--yukawa-type",
        type=int,
        default=None,
        help="Tipo de Yukawa (1..4, opcional).",
    )
    parser.add_argument(
        "--mphi-min",
        type=float,
        default=130.0,
        help="m_phi mínimo global (por defecto 130).",
    )
    parser.add_argument(
        "--mphi-max",
        type=float,
        default=290.0,
        help="m_phi máximo global (por defecto 290).",
    )
    parser.add_argument(
        "--n-chunks",
        type=int,
        default=16,
        help="Número de sub-intervalos (chunks) en m_phi (por defecto 16).",
    )
    parser.add_argument(
        "--time-available-hours",
        type=float,
        default=None,
        help="Horas de cómputo disponibles en esta sesión (aprox, opcional).",
    )
    parser.add_argument(
        "--avg-chunk-minutes",
        type=float,
        default=None,
        help="Override manual: minutos promedio por chunk (para estimar tiempos).",
    )
    parser.add_argument(
        "--executable",
        type=str,
        default="./PhysParamScan",
        help="Ruta al ejecutable PhysParamScan (por defecto ./PhysParamScan).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Solo mostrar el plan y los comandos sin ejecutarlos.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Número de threads OpenMP a usar en PhysParamScan (por defecto: todos disponibles).",
    )
    return parser.parse_args()


def make_tb_label(tan_beta: float) -> str:
    # Representación compacta sin .0 al final
    s = f"{tan_beta:g}"
    return s


def compute_chunks(mphi_min: float, mphi_max: float, n_chunks: int):
    if n_chunks <= 0:
        raise ValueError("n_chunks debe ser > 0")
    L = mphi_max - mphi_min
    if L <= 0:
        raise ValueError("Se requiere mphi_max > mphi_min")
    step = L / n_chunks
    chunks = []
    for i in range(n_chunks):
        left = mphi_min + i * step
        if i == n_chunks - 1:
            right = mphi_max
        else:
            right = mphi_min + (i + 1) * step
        chunks.append({"index": i, "mphi_min": left, "mphi_max": right})
    return chunks


def detect_completed_chunks(outdir: Path, n_chunks: int):
    done = set()
    for i in range(n_chunks):
        fname = outdir / f"chunk_{i:02d}.csv"
        if fname.exists():
            done.add(i)
    return done


def load_timings(outdir: Path) -> Dict[int, float]:
    path = outdir / "timings.json"
    if not path.exists():
        return {}
    try:
        with path.open() as f:
            data = json.load(f)
        return {int(k): float(v) for k, v in data.items()}
    except Exception:
        return {}


def save_timings(outdir: Path, timings: Dict[int, float]):
    path = outdir / "timings.json"
    serializable = {str(k): float(v) for k, v in timings.items()}
    with path.open("w") as f:
        json.dump(serializable, f, indent=2, sort_keys=True)


def estimate_avg_chunk_minutes(args, timings: Dict[int, float]):
    if args.avg_chunk_minutes is not None:
        return args.avg_chunk_minutes
    if timings:
        vals = list(timings.values())
        avg_sec = sum(vals) / len(vals)
        return avg_sec / 60.0
    return None


def main():
    args = parse_args()

    tb_label = make_tb_label(args.tan_beta)
    outdir = Path("outcsv") / f"tb{tb_label}"
    outdir.mkdir(parents=True, exist_ok=True)

    # Detectar threads disponibles
    available_threads = os.cpu_count()
    if available_threads is None:
        available_threads = 4
        print("[WARN] No se pudo detectar el número de threads disponibles. Usando valor por defecto: 4")

    # Determinar threads a usar
    selected_threads = args.threads if args.threads is not None else available_threads
    if selected_threads <= 0:
        print(f"[ERROR] --threads debe ser > 0. Usando {available_threads}.")
        selected_threads = available_threads

    print(f"[INFO] tan_beta = {args.tan_beta} (tb_label = tb{tb_label})")
    print(f"[INFO] Directorio de salida: {outdir}")
    print(f"[INFO] Threads disponibles: {available_threads}, a usar: {selected_threads}")

    # Definir chunks globales
    chunks = compute_chunks(args.mphi_min, args.mphi_max, args.n_chunks)

    # Detectar qué chunks ya tienen CSV
    done_chunks = detect_completed_chunks(outdir, args.n_chunks)
    pending_indices = [c["index"] for c in chunks if c["index"] not in done_chunks]
    pending_count = len(pending_indices)

    print(f"[INFO] Chunks totales: {args.n_chunks}")
    print(f"[INFO] Chunks completados: {len(done_chunks)} -> {sorted(done_chunks) if done_chunks else 'ninguno'}")
    print(f"[INFO] Chunks pendientes: {pending_count} -> {pending_indices if pending_indices else 'ninguno'}")

    if pending_count == 0:
        print("[INFO] Nada pendiente: el barrido ya está completo para este tan_beta.")
        return

    # Cargar tiempos históricos
    timings = load_timings(outdir)
    avg_chunk_min = estimate_avg_chunk_minutes(args, timings)

    suggested_chunks = pending_count
    est_total_min = None

    if args.time_available_hours is not None and avg_chunk_min is not None:
        capacity_min = args.time_available_hours * 60.0
        max_chunks_time = max(1, int(capacity_min // avg_chunk_min))
        suggested_chunks = max(1, min(max_chunks_time, pending_count))
        est_total_min = suggested_chunks * avg_chunk_min

        print(
            f"[INFO] Tiempo medio estimado por chunk: {avg_chunk_min:.1f} min "
            f"(basado en timings previos o override)."
        )
        print(
            f"[INFO] Con {args.time_available_hours:.2f} h se podrían ejecutar ~{max_chunks_time} chunks "
            f"(limitado a {suggested_chunks} por los pendientes)."
        )
        print(
            f"[INFO] Ejecutar {suggested_chunks} chunks tomaría aproximadamente "
            f"{est_total_min:.1f} minutos."
        )
    elif args.time_available_hours is not None and avg_chunk_min is None:
        print(
            "[WARN] Se indicó --time-available-hours pero no hay estimación de tiempo por chunk "
            "(ni override ni timings previos). No se puede estimar la duración total."
        )

    default_nrun = suggested_chunks if suggested_chunks > 0 else pending_count

    # Interacción para elegir cuántos chunks correr
    while True:
        try:
            raw = input(
                f"¿Cuántos chunks quieres ejecutar ahora? [1..{pending_count}] "
                f"(Enter = {default_nrun}): "
            ).strip()
        except EOFError:
            # En caso de no-interactivo, usar default
            raw = ""
        if not raw:
            nrun = default_nrun
        else:
            try:
                nrun = int(raw)
            except ValueError:
                print("Por favor ingresa un entero válido.")
                continue
        if 1 <= nrun <= pending_count:
            break
        print(f"Valor fuera de rango. Debe estar entre 1 y {pending_count}.")

    # Elegir los primeros nrun pendientes (estrategia: orden creciente)
    pending_indices.sort()
    to_run_indices = pending_indices[:nrun]

    # Recalcular estimación de tiempo para la selección final
    if avg_chunk_min is not None:
        est_total_min = nrun * avg_chunk_min
        print(
            f"[RESUMEN] Se planea ejecutar {nrun} chunks: {to_run_indices} "
            f"(tiempo estimado ~{est_total_min:.1f} min)."
        )
    else:
        print(f"[RESUMEN] Se planea ejecutar {nrun} chunks: {to_run_indices} (sin estimación de tiempo).")

    # Avisar si se completaría el barrido
    if len(done_chunks) + nrun == args.n_chunks:
        print("[INFO] Si todo va bien, con esta sesión se completará el barrido (todos los chunks).")
    else:
        print(
            f"[INFO] Tras esta sesión quedarán "
            f"{args.n_chunks - (len(done_chunks) + nrun)} chunks pendientes."
        )

    # Selección de threads (antes de confirmación final)
    while True:
        try:
            threads_input = input(
                f"¿Threads a usar? [1..{available_threads}] (Enter = {selected_threads}): "
            ).strip()
        except EOFError:
            # En caso de no-interactivo, usar selected_threads
            threads_input = ""
        
        if not threads_input:
            # Usar el valor actual
            break
        else:
            try:
                user_threads = int(threads_input)
                if 1 <= user_threads <= available_threads:
                    selected_threads = user_threads
                    print(f"[OK] Se usarán {selected_threads} threads.")
                    break
                else:
                    print(f"[WARN] Valor fuera de rango [1..{available_threads}]. Intenta de nuevo.")
            except ValueError:
                print("Por favor ingresa un entero válido.")

    # Confirmación final
    try:
        confirm = input("¿Confirmas la ejecución? [y/N]: ").strip().lower()
    except EOFError:
        confirm = "n"

    if confirm not in ("y", "yes", "s", "si", "sí"):
        print("[ABORTADO] No se ejecutó ningún comando.")
        return

    # Ejecutar los chunks seleccionados
    start_all = time.time()
    for idx in to_run_indices:
        ch = next(c for c in chunks if c["index"] == idx)
        mmin = ch["mphi_min"]
        mmax = ch["mphi_max"]
        outcsv = outdir / f"chunk_{idx:02d}.csv"

        cmd = [
            args.executable,
            str(outcsv),
            str(args.N_mphi),
            str(args.tan_beta),
        ]
        if args.tol_exp is not None:
            cmd.extend(["--tol_exp", str(args.tol_exp)])
        if args.N_m12 is not None:
            cmd.extend(["--N_m12", str(args.N_m12)])
        if args.yukawa_type is not None:
            cmd.extend(["--yukawa_type", str(args.yukawa_type)])

        cmd.extend(
            [
                "--mphi_min",
                f"{mmin:.6f}",
                "--mphi_max",
                f"{mmax:.6f}",
                "--N_mphi",
                "1",
                "--threads",
                str(selected_threads)
            ]
        )

        print("\n[RUN]", " ".join(cmd))
        t0 = time.time()
        if not args.dry_run:
            subprocess.run(cmd, check=True)
        t1 = time.time()
        elapsed = t1 - t0
        print(f"[DONE] chunk {idx:02d} en {elapsed/60.0:.2f} min")

        # Actualizar timings
        timings[idx] = elapsed
        save_timings(outdir, timings)

    total_elapsed = time.time() - start_all
    print(f"\n[FIN] Tiempo total de esta sesión: {total_elapsed/60.0:.1f} min")

    done_chunks = detect_completed_chunks(outdir, args.n_chunks)
    if len(done_chunks) == args.n_chunks:
        print("[FIN] Barrido COMPLETO para este tan_beta.")
    else:
        print(
            f"[FIN] Barrido PARCIAL. Chunks completados: {len(done_chunks)}/{args.n_chunks}. "
            f"Pendientes: {sorted(set(range(args.n_chunks)) - done_chunks)}"
        )


if __name__ == "__main__":
    main()
