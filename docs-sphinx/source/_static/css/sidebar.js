document.addEventListener("DOMContentLoaded", function () {
    const body = document.body;
    const button = document.createElement("button");

    button.className = "sidebar-toggle";
    button.type = "button";
    button.setAttribute("aria-label", "Toggle navigation sidebar");
    button.setAttribute("aria-expanded", "true");
    button.innerHTML = "☰";

    document.body.appendChild(button);

    /*
     * Restore the state saved on the previous page.
     * The sidebar is open by default when no state has been saved.
     */
    const sidebarCollapsed =
        localStorage.getItem("sidebar-collapsed") === "true";

    body.classList.toggle("sidebar-collapsed", sidebarCollapsed);
    button.setAttribute(
        "aria-expanded",
        String(!sidebarCollapsed)
    );

    button.addEventListener("click", function () {
        const collapsed =
            body.classList.toggle("sidebar-collapsed");

        localStorage.setItem(
            "sidebar-collapsed",
            String(collapsed)
        );

        button.setAttribute(
            "aria-expanded",
            String(!collapsed)
        );
    });
});
