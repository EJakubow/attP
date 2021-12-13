function startSearch() {
    const headsup = document.getElementById("headsup");
    headsup.style.display = "none";

    const start = document.getElementById("start");
    start.style.display = "none";

    const inProgress = document.getElementById("inProgress");
    inProgress.style.display = "block";

    window.location.href = "/in_progress";
}

$(document).ready(function() {
    $.each($('.navbar').find('li'), function() {
        const pathname = window.location.pathname;
        $(this).toggleClass('active',
            pathname.indexOf($(this).find('a').attr('href')) > -1);
        if(pathname === '/') {
            $('.navbar .search').addClass("active")
        };
    });
});
