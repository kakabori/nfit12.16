proc open_funclmdif_window {} {
  global flw
  toplevel .funclmdif
  set flw .funclmdif
  wm title $flw "funclmdif"
  text $flw.text -state normal -width 80 -height 10 -yscrollcommand "$flw.scroll set" -cursor arrow
  scrollbar $flw.scroll -command "$flw.text yview"
  pack $flw.scroll $flw.text -side right -fill y
}


proc show_warning {msg} {
  global flw
	$flw.text insert end $msg
	$flw.text see end
}
