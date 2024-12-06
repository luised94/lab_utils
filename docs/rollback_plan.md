ROLLBACK PLAN:
If issues occur at any step:
bash
# Record current commit for reference
git rev-parse HEAD > /tmp/last_working_commit

# Rollback
git reset --hard HEAD~1
./bash/tests/core/run_tests.sh

# If still broken
git reset --hard $(cat /tmp/last_working_commit)
