-- lw.lua
-- Neovim keybindings for lw (lab workflow)
-- Loaded by extension loader after lazy.nvim setup

-- === GUARDS ===
if vim == nil then
    print("lw.lua: not running in neovim, exiting")
    return
end

if vim.keymap == nil then
    vim.notify("lw.lua: vim.keymap not available (neovim < 0.7?), exiting", vim.log.levels.ERROR)
    return
end

if vim.api == nil then
    vim.notify("lw.lua: vim.api not available, exiting", vim.log.levels.ERROR)
    return
end

local has_telescope, telescope_builtin = pcall(require, "telescope.builtin")
if not has_telescope then
    vim.notify("[lw] telescope.builtin unavailable; install nvim-telescope/telescope.nvim", vim.log.levels.WARN)
    return nil
end

-- === CONFIGURATION ===
local api = vim.api
local fn  = vim.fn

local lw_lab_root = os.getenv("LW_LAB_ROOT")
if lw_lab_root == nil or lw_lab_root == "" then
    vim.notify("[lw] LW_LAB_ROOT not set. Source lw.bash first.", vim.log.levels.ERROR)
    return nil
end
lw_lab_root = fn.fnamemodify(lw_lab_root, ":p"):gsub("/$", "")

local lw_db    = os.getenv("LW_DB") or string.format("%s/lab.db", lw_lab_root)
local lw_notes = os.getenv("LW_NOTES") or string.format("%s/lab_notes.md", lw_lab_root)

local experiments_dir = string.format("%s/experiments", lw_lab_root)
local protocols_dir   = string.format("%s/protocols", lw_lab_root)

-- === CONSTANTS ===

---@type string
local MODE_NORMAL = "n"

---@type string
local LEADER_L = "<leader>l"

---@type string
local FIELD_SEP = "\t"

---@type string
local TODO_PATTERN = "TODO"

---@type string[]
local LAB_EXCLUDE_DIRS = {
    ".git",
    "__pycache__",
    ".venv",
    "staging",
}

-- === FUNCTIONS ===

---@param sql string
---@return string[]
local function query_db(sql)
    if fn.filereadable(lw_db) ~= 1 then
        vim.notify("[lw] database not found: " .. lw_db, vim.log.levels.ERROR)
        return {}
    end
    local cmd = string.format('sqlite3 -separator "%s" "%s" "%s"', FIELD_SEP, lw_db, sql)
    local raw = fn.system(cmd)
    if vim.v.shell_error ~= 0 then
        vim.notify("[lw] sqlite3 query failed", vim.log.levels.ERROR)
        return {}
    end
    return vim.split(raw, "\n", { trimempty = true })
end

---@param path string
---@return function
local function make_open_fn(path)
    return function()
        vim.cmd(string.format("edit %s", fn.fnameescape(path)))
    end
end

---@return nil
local function prepend_date_header()
    local bufnr = 0
    if not api.nvim_buf_get_option(bufnr, "modifiable") then
        return
    end

    local bufname = api.nvim_buf_get_name(bufnr)
    if bufname:match("lab_notes%.md$") == nil then
        vim.notify("[lw] date header only in lab_notes.md", vim.log.levels.WARN)
        return
    end

    local date_string = os.date("## %Y-%m-%d")
    local header_region_limit = 50
    local lines = api.nvim_buf_get_lines(bufnr, 0, header_region_limit, false)

    local existing_row = nil
    for i = 1, #lines do
        local line = lines[i]
        if line == date_string then
            existing_row = i
            break
        end
        if line:sub(1, 3) ~= "## " then
            break
        end
    end

    if existing_row ~= nil then
        vim.notify("[lw] today's header exists", vim.log.levels.INFO)
        api.nvim_win_set_cursor(0, { existing_row + 1, 0 })
        return
    end

    api.nvim_buf_set_lines(bufnr, 0, 0, false, { date_string, "" })
    api.nvim_win_set_cursor(0, { 2, 0 })
end

---@return nil
local function pick_experiment()
    local rows = query_db(
        "SELECT id, type, title, status FROM experiments ORDER BY date_started DESC"
    )
    if #rows == 0 then
        vim.notify("[lw] no experiments found", vim.log.levels.INFO)
        return
    end

    local entries = {}
    for _, row in ipairs(rows) do
        local parts = vim.split(row, FIELD_SEP, { plain = true })
        if #parts >= 4 then
            table.insert(entries, {
                id     = parts[1],
                etype  = parts[2],
                title  = parts[3],
                status = parts[4],
            })
        end
    end

    if #entries == 0 then
        vim.notify("[lw] failed to parse experiment rows", vim.log.levels.ERROR)
        return
    end

    ---@param item table
    ---@return string
    local function format_item(item)
        return string.format(
            "%s  %-14s %-24s [%s]",
            item.id, item.etype, item.title, item.status
        )
    end

    local select_opts = {
        prompt      = "Select experiment:",
        format_item = format_item,
    }

    ---@param choice table|nil
    ---@return nil
    local function on_selection(choice)
        if choice == nil then
            vim.notify("[lw] no experiment selected", vim.log.levels.INFO)
            return
        end
        api.nvim_put({ choice.id }, "c", true, true)
    end

    vim.ui.select(entries, select_opts, on_selection)
end

-- pick_strain (queried the retired `strains` table) parked in
-- archive/satellites_dormant.md; revive when strains rejoins LIVE_TABLE_NAMES.

---@return nil
local function find_lab_files()
    local file_ignore_patterns = {}
    for i = 1, #LAB_EXCLUDE_DIRS do
        local dir     = LAB_EXCLUDE_DIRS[i]
        local escaped = dir:gsub("([^%w])", "%%%1")
        file_ignore_patterns[#file_ignore_patterns + 1] = string.format("^%s/", escaped)
    end

    telescope_builtin.find_files({
        prompt_title         = string.format("Lab (%s)", lw_lab_root),
        cwd                  = lw_lab_root,
        hidden               = false,
        follow               = true,
        file_ignore_patterns = file_ignore_patterns,
    })
end

---@return nil
local function find_protocols()
    telescope_builtin.find_files({
        prompt_title = "Protocols",
        cwd          = protocols_dir,
        hidden       = false,
        follow       = true,
    })
end

---@return nil
local function find_experiment_files()
    local bufname = api.nvim_buf_get_name(0)
    local exp_match = bufname:match("experiments/(%d%d%d%d%d%d%d%d_LM%-%d%d%d%d_[^/]+)/")

    if exp_match == nil then
        -- not inside an experiment folder; prompt for experiment
        local rows = query_db(
            "SELECT id, folder_name, title FROM experiments WHERE status = 'active' "
            .. "ORDER BY date_started DESC"
        )
        if #rows == 0 then
            vim.notify("[lw] no active experiments", vim.log.levels.INFO)
            return
        end

        local entries = {}
        for _, row in ipairs(rows) do
            local parts = vim.split(row, FIELD_SEP, { plain = true })
            if #parts >= 3 then
                table.insert(entries, {
                    id     = parts[1],
                    folder = parts[2],
                    title  = parts[3],
                })
            end
        end

        local function format_item(item)
            return string.format("%s  %s", item.id, item.title)
        end

        vim.ui.select(entries, {
            prompt      = "Select experiment folder:",
            format_item = format_item,
        }, function(choice)
            if choice == nil then return end
            local folder_path = string.format("%s/%s", experiments_dir, choice.folder)
            telescope_builtin.find_files({
                prompt_title = string.format("%s (%s)", choice.id, choice.title),
                cwd          = folder_path,
                hidden       = false,
                follow       = true,
            })
        end)
        return
    end

    -- inside an experiment folder; scope telescope to it
    local folder_path = string.format("%s/%s", experiments_dir, exp_match)
    local exp_id = exp_match:match("(LM%-%d%d%d%d)")
    telescope_builtin.find_files({
        prompt_title = string.format("%s files", exp_id or exp_match),
        cwd          = folder_path,
        hidden       = false,
        follow       = true,
    })
end

---@return nil
local function lab_todos()
    local files_to_search = { lw_notes }
    local todo_entries     = {}

    for _, file_path in ipairs(files_to_search) do
        if fn.filereadable(file_path) ~= 1 then
            goto continue_file
        end
        do
            local file_lines     = fn.readfile(file_path)
            local current_header = "(no section)"
            for line_number, line_text in ipairs(file_lines) do
                if line_text:match("^## %d%d%d%d%-%d%d%-%d%d") then
                    current_header = line_text
                end
                if line_text:match(TODO_PATTERN) then
                    table.insert(todo_entries, {
                        filename = file_path,
                        lnum     = line_number,
                        text     = string.format("%s | %s", current_header, vim.trim(line_text)),
                    })
                end
            end
        end
        ::continue_file::
    end

    if #todo_entries == 0 then
        vim.notify("[lw] no TODOs found", vim.log.levels.INFO)
        return
    end

    fn.setqflist({}, " ", {
        title = "Lab TODOs",
        items = todo_entries,
    })
    vim.cmd("copen")
    vim.notify(
        string.format("[lw] %d TODOs loaded into quickfix", #todo_entries),
        vim.log.levels.INFO
    )
end

-- === DECLARATIONS ===

---@type table[]
local commands = {
    {
        name = "LwPickExperiment",
        fn   = pick_experiment,
        opts = { desc = "lw: pick experiment and insert ID" },
    },
    {
        name = "LwTodos",
        fn   = lab_todos,
        opts = { desc = "lw: list TODOs in quickfix" },
    },
}

---@class LwKeymapSpec
---@field key  string
---@field fn   function
---@field desc string

---@type LwKeymapSpec[]
local LW_KEYMAP_SPECS = {
    -- file openers
    { key = "n", fn = make_open_fn(lw_notes), desc = "open lab notes"                   },
    -- actions
    { key = "h", fn = prepend_date_header,    desc = "insert date header in lab notes"   },
    { key = "e", fn = pick_experiment,        desc = "pick experiment, insert ID"        },
    -- 's' (strain pick) parked in archive/satellites_dormant.md; strains table retired.
    -- navigation
    { key = "f", fn = find_lab_files,         desc = "find files in lab root (telescope)" },
    { key = "p", fn = find_protocols,         desc = "find protocols (telescope)"        },
    { key = "x", fn = find_experiment_files,  desc = "find files in experiment folder"   },
    -- review
    { key = "t", fn = lab_todos,              desc = "list TODOs in quickfix"            },
}

local keymaps   = {}
local seen_keys = {}
for _, spec in ipairs(LW_KEYMAP_SPECS) do
    if spec.fn == nil then
        vim.notify(string.format("[lw] keymap '%s': nil function, skipping", spec.key), vim.log.levels.WARN)
    elseif seen_keys[spec.key] then
        vim.notify(string.format("[lw] keymap '%s': duplicate key, skipping", spec.key), vim.log.levels.WARN)
    else
        seen_keys[spec.key] = true
        keymaps[#keymaps + 1] = { MODE_NORMAL, LEADER_L .. spec.key, spec.fn, { desc = "lw: " .. spec.desc } }
    end
end

return {
    keymaps  = keymaps,
    autocmds = {},
    commands = commands,
    setup    = nil,
}
