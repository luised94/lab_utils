# Dormant satellite code

Parked code pulled out of the live `lw` satellite files (`lw.sh`, `lw.lua`)
when the rebuild collapsed the tool to its current 6-command / 2-table
surface. None of this is dead in the archaeological sense: each section names
the condition under which it returns. Revive a section by copy-pasting it back
into its source file, then delete the section here.

Source of truth for the live surface is `lw.py` (`COMMAND_USAGE` for commands,
`LIVE_TABLE_NAMES` for tables).

---

## lw.sh completion: pre-migration subcommand tree

Parked 2026-06-26. Revive when `exp` / `strain` / `stage` subcommands return to
`lw.py`. Verify still-absent: `rg 'COMMAND_USAGE' lw.py` lists the live
commands; if it still shows only init/new/link/show/list/complete, this tree
stays parked.

The original `_lw_completions` advertised a command tree that no longer exists
(tab-completing to commands that return usage errors is worse than no
completion). Old sub-tables:

```bash
    commands="init status exp strain stage"
    exp_subs="init show update complete link addstrain delete find list manifest"
    strain_subs="add show update list"
    stage_subs="list assign expect auto"

    if [ "${COMP_CWORD}" -eq 1 ]; then
        COMPREPLY=($(compgen -W "${commands}" -- "${cur}"))
    elif [ "${COMP_CWORD}" -eq 2 ]; then
        case "${prev}" in
            exp)    COMPREPLY=($(compgen -W "${exp_subs}" -- "${cur}")) ;;
            strain) COMPREPLY=($(compgen -W "${strain_subs}" -- "${cur}")) ;;
            stage)  COMPREPLY=($(compgen -W "${stage_subs}" -- "${cur}")) ;;
        esac
    fi
```

---

## lw.lua strain picker

Parked 2026-06-26. Queries the retired `strains` table. Revive when `strains`
rejoins `LIVE_TABLE_NAMES` in `lw.py`. Verify still-dead: `rg 'strains' lw.py`
should appear only in the `schema_notes` orphaned-tables list.

When reviving: paste `pick_strain` back among the functions, re-add the
`LwPickStrain` entry to the `commands` table, and re-add the
`{ key = "s", fn = pick_strain, ... }` entry to `LW_KEYMAP_SPECS`.

```lua
---@return nil
local function pick_strain()
    local rows = query_db(
        "SELECT id, COALESCE(label, '*NO LABEL*'), genotype FROM strains ORDER BY id"
    )
    if #rows == 0 then
        vim.notify("[lw] no strains found", vim.log.levels.INFO)
        return
    end

    local entries = {}
    for _, row in ipairs(rows) do
        local parts = vim.split(row, FIELD_SEP, { plain = true })
        if #parts >= 3 then
            table.insert(entries, {
                id       = parts[1],
                label    = parts[2],
                genotype = parts[3],
            })
        end
    end

    if #entries == 0 then
        vim.notify("[lw] failed to parse strain rows", vim.log.levels.ERROR)
        return
    end

    ---@param item table
    ---@return string
    local function format_item(item)
        local geno = item.genotype
        if #geno > 40 then
            geno = geno:sub(1, 37) .. "..."
        end
        return string.format("%-12s %-16s %s", item.id, item.label, geno)
    end

    local select_opts = {
        prompt      = "Select strain:",
        format_item = format_item,
    }

    ---@param choice table|nil
    ---@return nil
    local function on_selection(choice)
        if choice == nil then
            vim.notify("[lw] no strain selected", vim.log.levels.INFO)
            return
        end
        api.nvim_put({ choice.id }, "c", true, true)
    end

    vim.ui.select(entries, select_opts, on_selection)
end
```

Command-table entry (was in `commands`):

```lua
    {
        name = "LwPickStrain",
        fn   = pick_strain,
        opts = { desc = "lw: pick strain and insert ID" },
    },
```

Keymap entry (was in `LW_KEYMAP_SPECS`):

```lua
    { key = "s", fn = pick_strain,            desc = "pick strain, insert ID"            },
```
