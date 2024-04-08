# 1. Organize data in a single directory with sub-directories corresponding to each data type (ATAC, RNA, etc.)
# Important: GEO said to NOT upload metadata through FTP!!

# 2. When data is ready to upload, log onto to GEO: https://www.ncbi.nlm.nih.gov/geo/submitter/

# 3. Go to "GEO submissions" or "My submissions" in upper-right taskbar

# 4. Click "New submission"

# 5. Click "Submit high-throughput sequencing"

# 6. Under the "Uploading your submission" section, click "Transfer Files".

# 7. On the terminal in flh1, open a screen session.
screen

# 8. Upload data onto GEO using lftp. You will need the following info:
# the GEO-provided password under "Step 2", bullet point f
# the GEO-provided upload space under "Step 1"
# The path to your submission directory

lftp ftp://geoftp:[password]>@ftp-private.ncbi.nlm.nih.gov
cd [upload space]
mirror -R [path_to_submission_directory]

# 9. Exit screen using Ctrl+A+D

# 10. To go back to your screen session to check progress
# If you forgot your screen_session_id, just run screen R and it'll list out all sessions that are active
screen -R [screen_session_id]
