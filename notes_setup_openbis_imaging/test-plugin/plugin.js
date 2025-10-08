function CustomPluginViewTechnology() {
    this.init();
}

$.extend(CustomPluginViewTechnology.prototype, ELNLIMSPlugin.prototype, {

    init: function() {

    },

    getExtraUtilities : function() {
        var _this = this;
        return [{
                icon : "fa fa-info-circle",
                uniqueViewName : "TEMPLATE_EXTRA_UTILITIES_CUSTOM_VIEW_JS_TEMPLATE",
                label : "AiiDAlab",
                paintView : function($header, $content) {
                $header.append($("<h1>").append("AiiDAlab Link"));
                                    var src = "http://localhost:8888/apps/apps/home/start.ipynb";
                                    var win = window.open(src, '_blank');
                                    win.focus();
                }
            }
        ];
    }
});

profile.plugins.push(new CustomPluginViewTechnology());
