"""
prepare_outputfolder(outputfolder)

    This function checks if the outputfolder exists. If it doesn't, it creates it.
    If it does exist, it deletes the folder and its contents, then recreates it.
    Please note that this function will delete all files and subdirectories in the
    specified outputfolder, so use it with caution.
"""
"""
    prepare_outputfolder(outputfolder)

Make `outputfolder` ready to be written into, moving any existing one aside to
`<name>_previous` rather than deleting it. Starting from an empty folder is what keeps
stale plots and data from a previous fit out of the current results; keeping the old one is
what stops a mistyped folder name destroying its contents.
"""
prepare_outputfolder(outputfolder) = backupfolder(outputfolder)