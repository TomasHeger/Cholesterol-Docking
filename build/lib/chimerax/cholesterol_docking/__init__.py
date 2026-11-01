from chimerax.core.toolshed import BundleAPI

class _CholesterolDockingAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bi, ci, logger):
        if ci.name != "choldock7":
            raise ValueError("Unknown command: %s" % ci.name)
        from .cmd import choldock7, choldock7_desc
        if choldock7_desc.synopsis is None:
            choldock7_desc.synopsis = ci.synopsis
        from chimerax.core.commands import register
        register(ci.name, choldock7_desc, choldock7, logger=logger)

    @staticmethod
    def start_tool(session, bi, ti):
        if ti.name != "Cholesterol Docking V7":
            raise ValueError("Unknown tool: %s" % ti.name)
        from .tool import CholesterolDockingTool
        return CholesterolDockingTool(session, ti.name)

bundle_api = _CholesterolDockingAPI()
