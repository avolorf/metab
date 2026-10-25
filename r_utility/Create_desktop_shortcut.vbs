Set sh = CreateObject("WScript.Shell")
d = CreateObject("Scripting.FileSystemObject").GetParentFolderName(WScript.ScriptFullName)
Set lnk = sh.CreateShortcut(sh.SpecialFolders("Desktop") & "\UroMetabR.lnk")
lnk.TargetPath = d & "\Start_UroMetabR.bat"
lnk.WorkingDirectory = d
lnk.Description = "UroMetabR"
lnk.Save
MsgBox "Shortcut UroMetabR created on Desktop"
