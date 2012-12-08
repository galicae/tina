void initStatus(Portallayout param) {
	List panelchildren = param.getChildren();
	for (int i = 0; i < panelchildren.size(); i++) {
		List panelIds = session.getAttribute("PortalChildren" + i);
		if (panelIds != null) {
			for (String panelId : panelIds) {
				Panel newPanel = param.getFellow(panelId);
				if (panelchildren.size() > 0)
					panelchildren.get(i).insertBefore(newPanel, panelchildren.get(0));
				else
					newPanel.setParent(panelchildren.get(i));
			}
		}
	}
}

/* Save Portlet Position (Session Scope) */
void saveStatus(Portallayout param) {
	int i = 0;
	for (Portalchildren ptlc : param.getChildren()) {
		List portletIds = new ArrayList();
		for (Panel portlet : ptlc.getChildren())
		portletIds.add(portlet.getId());
		session.setAttribute("PortalChildren" + i++, portletIds);
	}
}