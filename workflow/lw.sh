# ============================================================
# lw.sh -- lab workflow environment configuration
# ============================================================
# Source this file from bashrc or symlink into your sourced
# directory. All variables use the LW_ prefix.
#
# Derives LW_WINDOWS_USER from MC_WINDOWS_USER if set.
# Override any variable by exporting it before sourcing.
# ============================================================

# --- user configuration ---
LW_WINDOWS_USER="${LW_WINDOWS_USER:-${MC_WINDOWS_USER:-}}"
if [ -z "${LW_WINDOWS_USER}" ]; then
    echo "lw: error: Set MC_WINDOWS_USER or LW_WINDOWS_USER before sourcing lw.bash"
    return 1 2>/dev/null || true
fi

# --- editor (EDITOR -> nvim -> vim -> vi) ---
if [ -n "${EDITOR:-}" ] && command -v "${EDITOR}" >/dev/null 2>&1; then
    LW_EDITOR="${EDITOR}"
elif command -v nvim >/dev/null 2>&1; then
    LW_EDITOR="nvim"
elif command -v vim >/dev/null 2>&1; then
    LW_EDITOR="vim"
elif command -v vi >/dev/null 2>&1; then
    LW_EDITOR="vi"
else
    echo "lw: warning: No editor found (tried EDITOR, nvim, vim, vi). lwnotes will not work."
    LW_EDITOR=""
fi

# --- paths ---
LW_LAB_ROOT="/mnt/c/Users/${LW_WINDOWS_USER}/Desktop/lab"
LW_SCRIPT="/home/${USER}/personal_repos/lab_utils/workflow/lw.py"
LW_NOTES="${LW_LAB_ROOT}/lab_notes.md"
LW_DB="${LW_LAB_ROOT}/lab.db"

# --- backup paths ---
LW_DROPBOX_PATH="${LW_DROPBOX_PATH:-${MC_DROPBOX_PATH:-}}"
if [ -n "${LW_DROPBOX_PATH}" ]; then
    LW_BACKUP_DIR="${LW_DROPBOX_PATH}/lab-backup"
    LW_BACKUP_SCRIPT="/home/${USER}/personal_repos/lab_utils/workflow/lw_backup.sh"
else
    LW_BACKUP_DIR=""
    LW_BACKUP_SCRIPT=""
    echo "lw: warning: No Dropbox path set (MC_DROPBOX_PATH or LW_DROPBOX_PATH). lwbackup will not work."
fi

export LW_WINDOWS_USER
export LW_LAB_ROOT
export LW_SCRIPT
export LW_NOTES
export LW_DB
export LW_EDITOR
export LW_DROPBOX_PATH
export LW_BACKUP_DIR
export LW_BACKUP_SCRIPT

# --- dependency checks ---
if [ ! -f "${LW_SCRIPT}" ]; then
    echo "lw: warning: Script not found: ${LW_SCRIPT}"
fi

if [ ! -x "$HOME/.local/bin/uv" ]; then
    echo "lw: warning: uv not found at $HOME/.local/bin/uv. lw alias will not work."
fi

if [ ! -d "${LW_LAB_ROOT}" ]; then
    echo "lw: warning: Lab root does not exist: ${LW_LAB_ROOT}"
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
    echo "lw: warning: sqlite3 not found. lwcheck and lwbackup will not work."
fi

# --- aliases ---
alias lw='uv run "${LW_SCRIPT}"'

if [ -n "${LW_EDITOR}" ]; then
    alias lwnotes='${LW_EDITOR} "${LW_NOTES}"'
fi

# open the database for ad hoc queries
alias lwdb='sqlite3 "${LW_DB}"'

# database integrity check
alias lwcheck='sqlite3 "${LW_DB}" "PRAGMA integrity_check; PRAGMA foreign_key_check;"'

# backup to dropbox
if [ -n "${LW_BACKUP_SCRIPT}" ] && [ -f "${LW_BACKUP_SCRIPT}" ]; then
    alias lwbackup='bash "${LW_BACKUP_SCRIPT}"'
fi

# cd into lab root, or into an experiment folder by ID
lwcd() {
    if [ -z "$1" ]; then
        cd "${LW_LAB_ROOT}" || return 1
        return 0
    fi
    local match
    match=$(find "${LW_LAB_ROOT}/experiments" -maxdepth 1 -type d -name "*$1*" | head -1)
    if [ -n "${match}" ]; then
        cd "${match}" || return 1
    else
        echo "No experiment folder matching '$1'"
        return 1
    fi
}

# --- tab completion ---
# Mirrors lw.py's command surface. Verify: rg 'COMMAND_USAGE' lw.py.
# Update when commands are added or removed there. The pre-migration
# exp/strain/stage subcommand tree is parked in archive/satellites_dormant.md.
_lw_completions() {
    local cur prev commands
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    commands="init new link show list complete"

    if [ "${COMP_CWORD}" -eq 1 ]; then
        COMPREPLY=($(compgen -W "${commands}" -- "${cur}"))
    elif [ "${COMP_CWORD}" -eq 2 ]; then
        case "${prev}" in
            new)      COMPREPLY=($(compgen -W "--title" -- "${cur}")) ;;
            show)     COMPREPLY=($(compgen -W "--files" -- "${cur}")) ;;
            list)     COMPREPLY=($(compgen -W "--type --status" -- "${cur}")) ;;
        esac
    fi
}
complete -F _lw_completions lw
