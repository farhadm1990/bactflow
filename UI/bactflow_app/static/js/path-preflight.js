/* Shared directory preflight for QC / Assembly / Post-assembly. */
(function (global) {
  function formFrom(elOrForm) {
    if (!elOrForm) {
      return null;
    }
    if (elOrForm instanceof FormData) {
      return elOrForm;
    }
    if (typeof elOrForm === "string") {
      const node = document.getElementById(elOrForm);
      return node ? new FormData(node) : null;
    }
    if (elOrForm.tagName === "FORM") {
      return new FormData(elOrForm);
    }
    return null;
  }

  async function validatePaths(module, formLike) {
    const formData = formFrom(formLike);
    if (!formData) {
      return { ok: false, message: "Form not found for path validation." };
    }
    formData.set("module", module || "");
    let response;
    try {
      response = await fetch("/validate-paths", { method: "POST", body: formData });
    } catch (err) {
      return { ok: true, skipped: true, message: "" };
    }
    const text = ((await response.text()) || "").trim();
    let data = null;
    if (text) {
      try {
        data = JSON.parse(text);
      } catch (err) {
        data = null;
      }
    }
    if (!data && (response.status === 404 || /<!DOCTYPE|<html/i.test(text))) {
      return { ok: true, skipped: true, message: "" };
    }
    if (!response.ok || !data || !data.ok) {
      return {
        ok: false,
        message: (data && data.message) || text || "Directory path check failed.",
        fields: (data && data.fields) || [],
      };
    }
    return data;
  }

  function bindClearOnEdit(ids, onClear) {
    (ids || []).forEach((id) => {
      const el = document.getElementById(id);
      if (!el) {
        return;
      }
      const handler = () => {
        if (typeof onClear === "function") {
          onClear();
        }
      };
      el.addEventListener("input", handler);
      el.addEventListener("change", handler);
    });
  }

  global.BactflowPathPreflight = {
    validatePaths,
    bindClearOnEdit,
  };
})(window);
