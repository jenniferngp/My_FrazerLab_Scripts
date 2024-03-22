
# 0. Make the group repository public (already done!)

# 1. Fork repo into your private account

# 2. Generate a fine-grained token that allows you to commit remotely
# - Click on your profile icon (top right)
# - Go to "Settings"
# - Select "Developer Settings" on the left navigation bar
# - Select "Personal access tokens"
# - Select "Fine-grained tokens"
# - Select "Generate new token". Complete the form. Under "Repository access", I selected "Only select repositories". Under "Repository permissions", make sure you enable "Read and Write" access for "Actions", "Administration", "Commit statuses", "Contents", "Pull requests", "Deployments" (I think optional). 
# - Once your token is generated, save the password to a different location. You will not see it again

# 3. On the terminal, git clone the forked repo. I did this in my own personal directory. 
git clone https://github.com/<username>/iPSCORE_QTL_Resource.git

# 4. Add/delete/update any files on the terminal, in that newly cloned directory.

# 5. To upload changes, create a new branch
git checkout -b new-branch # creates a new branch "new-branch"
git branch # checks which branch you're in

# 6. Stage changes
# https://stackoverflow.com/questions/41497031/add-all-removed-files-to-a-commit-with-git
git status # lists what changes have been made
git add -A # "-A" stages all changes (new files, modified, and deleted)
git status # check that all changes have been staged and ready to commit

# 7. Commit with a detailed message
git commit -m "message" 
git status # check that there is "nothing to commit, working directory clean"

# 9. Push changes 
git push origin new-branch
# Then, enter your github username and the token password you generated in step 2. 

# 10. On the brower, go to the forked repo in your account. 

# 11. Select "Compare & pull request", then "Create pull request".

# 12. Since you are already the manager of the original repo, 
# the page will take you to the original repo page and ask if you want to merge the pull request. 
# If yes, select "Merge pull request", then "Confirm merge", then "Delete branch" (unless you have more to commit). 

# 13. IMPORTANT! Sync the repo on your terminal with the main branch! 
# IMPORTANT!! Do this step AGAIN before you commit next time in case the repo has been updated by someone else!
git pull origin main

# 14. IMPORTANT!! On the browser, go to your forked repo on GitHub and select "Sync fork"! 
# IMPORTANT!! Do this step AGAIN before you commit next time in case the repo has been updated by someone else!

# 15. Done! 

